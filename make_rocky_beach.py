import xarray as xr
import rioxarray as xrr
import logging
from shared import *

def get_enso_poly(gdf_enso_map_in, da_dem_in, poly_in):
    """ Get the ENSO dem values within a polygon """
    gdf_enso_map_poly = gdf_enso_map_in[gdf_enso_map_in.contains(poly_in)]
    if gdf_enso_map_poly.empty:
        return 0
        
    elif gdf_enso_map_poly.shape[0]>1:
        return gdf_enso_map_poly.shape[0]

    path = gdf_enso_map_poly['path'].item()
    da_enso_dem = xrr.open_rasterio(path)
    da_enso_dem_re = da_enso_dem.rio.reproject(da_dem_in.rio.crs).sel(band=1)
    da_enso_dem_re = da_enso_dem_re.where(da_enso_dem_re>-1e20)
    da_enso_dem_re.rio.write_nodata(da_enso_dem_re.rio.nodata, encoded=True, inplace=True)
    # downsample from 0.5 m to 1 m horizontal resolution to match CoNED and MLLW/MAH
    da_enso_dem_re1 = da_enso_dem_re.interp_like(da_dem_in, method='linear')
    # clip to polygon
    try:
        da_enso_poly = da_enso_dem_re1.rio.clip([poly_in], da_dem_in.rio.crs, all_touched=True)
    except Exception as E:
        log.error(E)
        da_enso_poly = path
    return da_enso_poly


def compare_elevations_poly(poly_in, da_dem_in, da_mw_0_in, da_mw_slr_in, kind='MLLW'):
    """ The actual work: compare the DEM elevations within the poly against  current/future SLR """
    kind = kind.upper()
    assert kind in 'MLLW MAH MAH_ENSO'.split(), f'Incorrect kind: {kind}'
    scen = da_mw_slr_in.attrs['scenario']

    # get the mw in this polygon
    da_mw_0_poly = da_mw_0_in.rio.clip([poly_in], da_mw_0_in.rio.crs, all_touched=True)
    da_mw_slr_poly = da_mw_slr_in.rio.clip([poly_in], da_mw_slr_in.rio.crs, all_touched=True)
    
    if kind == 'MLLW':
        da_dem_mw_0_poly = (da_dem_in>da_mw_0_poly).astype(int)
        da_dem_mw_slr_poly = (da_dem_in>da_mw_slr_poly).astype(int)
 
        # preserve nodata values in dem for more accurate percentages
        da_dem_mw_0_poly = da_dem_mw_0_poly.where(~da_dem_in.isnull())
        da_dem_mw_slr_poly = da_dem_mw_slr_poly.where(~da_dem_in.isnull())
        
        pct0  = 100*(da_dem_mw_0_poly.mean())
        pct1  = 100*(da_dem_mw_slr_poly.mean())
        pct_lost = pct0-pct1
        msg = f'\t{pct0.item():.3f}% ABOVE {kind} today; {pct1.item():.3f} in {scen}. \tLost: {pct_lost.item():.3f}% ' 
        
    if kind in 'MAH MAH_ENSO'.split():
        # compare elevations against MAH --- ELEVATION BELOW MAH
        da_dem_mw_0_poly = (da_dem_in<da_mw_0_in).astype(int)
        da_dem_mw_slr_poly = (da_dem_in<da_mw_slr_in).astype(int)
        
        # preserve nodata values in dem for more accurate percentages
        da_dem_mw_0_poly = da_dem_mw_0_poly.where(~da_dem_in.isnull())
        da_dem_mw_slr_poly = da_dem_mw_slr_poly.where(~da_dem_in.isnull())
        pct0_mah  = 100*(da_dem_mw_0_poly.mean())
        pct1_mah  = 100*(da_dem_mw_slr_poly.mean())
        if pct0_mah > pct1_mah:
            pct_chg_mah = pct0_mah-pct1_mah
            msg = f'\t{pct0_mah.item():.3f}% BELOW {kind} today; {pct1_mah.item():.3f} in {scen}. \tLost: {pct_chg_mah.item():.3f}% ' 
        else:
            pct_chg_mah = pct1_mah-pct0_mah
            msg = f'\t{pct0_mah.item():.3f}% ABOVE {kind} today; {pct1_mah.item():.3f} in {scen}. \tGained: {pct_chg_mah.item():.3f}% ' 
            
    log.critical(msg)
    
    df_mw = da_dem_mw_0_poly.rename(kind).to_dataframe()[kind].reset_index().dropna()
    df_mw_slr = da_dem_mw_slr_poly.rename(scen).to_dataframe()[scen].reset_index().dropna()
    df_mw[f'{scen}_{kind}'] = df_mw_slr[scen]
            
    return df_mw


def main(region, habit='rocky', scen='med_rsl2050', path_wd=None, use_s2=False, test=False):
    habit = habit.lower()
    tstl = '_test' if test else ''
    s2 = '_s2' if use_s2 else ''
    assert habit in 'beaches rocky'.split(), f'Incorrect habit: {habit}'
    if habit == 'rocky':
        assert not use_s2, f'Sentinel2 not supported yet for {habit}'

    Obj  = SetupProj(region, habit, scen, path_wd, use_s2)
    path_log = Obj.path_res / f'log_{Obj.region.title()}_{habit}_pct{tstl}{s2}.txt'
    log2file(path_log)
    st0  = time.time()
    gdf_enso_map = get_enso_map(Obj.path_enso_dems, Obj.epsg) # enso dem path and its bounds
    
    lst_gsers = []
    # iterate over all DEM tiles
    for i, dem in enumerate((Obj.path_ca_dems).glob(f'{Obj.stem}.tif')):
        path_res = dem.parent / f'{dem.stem}_{habit}{tstl}{s2}.GeoJSON'
        # dem = Obj.path_ca_dems / 'CentCA_south_Topobathy_CoNED_1m_B2.tif' for testing
        if path_res.exists() and not test:
            log.info (f'gdf of {habit} in tile {dem.stem} exists, skipping...')
            continue
        
        sti = time.time()
        # get the tile 
        da_dem = xrr.open_rasterio(dem).sel(band=1)
        wi, si, ei, ni = da_dem.rio.bounds()

        da_dem = da_dem.where(da_dem>-3e4)
        da_dem.rio.write_nodata(da_dem.rio.nodata, encoded=True, inplace=True)

        # crop cari polygons within this tile
        gdf_cari_tile = Obj.gdf_cari0.cx[wi:ei, si:ni]
        n_polys = gdf_cari_tile.shape[0]
        if n_polys == 0:
            log.warning(f'No {habit} polygons within {dem.stem}, skipping tile.')
            continue
        else:
            log.info(f'{n_polys} {habit} polygons within {dem.stem}')
        
        # check there is actual coned data within the polygons
        da_dem_check = da_dem.rio.clip(gdf_cari_tile.geometry.tolist(), gdf_cari_tile.crs, 
                                       drop=True, invert=False, all_touched=True)
        if da_dem_check.isnull().all():
            log.warning(f'No actual DEM values within ANY {habit} polygons. Skipping coned tile: {dem.stem}.')
            continue
            
        # interpolate the mllw and mah to the 1m coned tiles
        da_mllw_re = Obj.da_mllw_0.interp_like(da_dem, method='linear')
        da_mllw_slr_re = Obj.da_mllw_slr.interp_like(da_dem, method='linear')

        da_mah_re = Obj.da_mah_0.interp_like(da_dem, method='linear')
        da_mah_slr_re = Obj.da_mah_slr.interp_like(da_dem, method='linear')
        
        elapi = time.time()-sti
        log.debug(f'Preparing DEM, MLLW, MAH took {elapi/60:.2f} min.')

        # iterate over each polygon in this tile
            ## this could probably be vectorized, as in da_dem_check
        lst_gdfs = []
        for j, (poly_ix, row) in enumerate(gdf_cari_tile.iterrows()):
            stj = time.time()
            poly = row.geometry
            cari_id = row['CARI_id']
            
            da_dem1 = da_dem.rio.clip([poly], gdf_cari_tile.crs, drop=True, invert=False, all_touched=True)
            if da_dem1.isnull().all():
                log.warning(f'No actual DEM values within {habit} polygon {poly_ix} of {dem.stem}, skipping poly_ix {poly_ix}.')
            
            print (''); log.critical (f'{dem.stem}, polyix={poly_ix} cari id={cari_id}:')
            df_mllw = compare_elevations_poly(poly, da_dem1, da_mllw_re, da_mllw_slr_re, 'MLLW')
            df_mah = compare_elevations_poly(poly, da_dem1, da_mah_re, da_mah_slr_re, 'MAH')
            
            da_enso_dem = get_enso_poly(gdf_enso_map, da_dem1, poly)
            if isinstance(da_enso_dem, int):
                log.warning (f'{da_enso_dem} ENSO dems with polygon cari_id={cari_id}, poly_ix={poly_ix}, skipping.') 
                continue
            elif isinstance(da_enso_dem, str): # actual returns the path 
                log.error(f'poly_ix={poly_ix}, cari_id={cari_id} not in ENSO DEM: {da_enso_dem}, skipping')
                continue
                
            df_mah_enso = compare_elevations_poly(poly, da_enso_dem, da_mah_re, da_mah_slr_re, 'MAH_ENSO')
            
            df_m = pd.merge(df_mllw, df_mah, on='x y'.split(), how='outer')
            df_m = pd.merge(df_m, df_mah_enso, on='x y'.split(), how='outer')
            df_m['poly_ix'] = poly_ix
            df_m['cari_id'] = cari_id
            
            ## enso may be smaller...
            # assert df_m.shape[0] == df_mllw.shape[0] == df_mah.shape[0] == df_mah_enso.shape[0], 'Shapes arent the same?'
            
            lst_gdfs.append(df2gdf(df_m, 'all', epsg=Obj.epsg))

            if test and j == 0:
                return da_dem1, da_mllw_re, da_mllw_slr_re, df_m, poly

            del df_m, da_dem1, df_mllw, df_mah, da_enso_dem, df_mah_enso
                
        gdf_tile = pd.concat(lst_gdfs)
        gdf_tile.to_file(path_res)
        log.critical(f'CONED tile completed: Wrote: {path_res}')
        del gdf_tile, lst_gdfs, da_dem, da_dem_check, gdf_cari_tile

    elap_tot = time.time()-st0
    log.critical (f'Finished writing epsg={Obj.epsg} CoNED {habit} Percents in {elap_tot/60:.2f} min.')
    log.info (f'Log file at:', path_log)
    return 


def concat_beaches(region, scen='med_rsl2050', path_wd=None, use_s2=False):
    """ Concatenate the beaches into a single (big) geojson file (e.g. 'South_beaches.GeoJSON') """
    Obj = CARIRegions(region, path_wd, scen, use_s2)
    s2_ext = '_s2' if use_s2 else ''
    paths_res = sorted(list(Obj.path_ca_dems.glob(f'*beach{s2_ext}.GeoJSON')))
    dst = Obj.path_res / f'{region}_beaches{s2_ext}.GeoJSON' 

    log.info (f'Concatenating: {len(paths_res)} {region} CoNED tiles')
    lst_gsers = []
    for i, path in enumerate(paths_res):
        log.info (f'Tile {i} of {len(paths_res)}') if i % 5 == 0  else ''
        gdfi  = gpd.read_file(path)
        gdfi[scen.replace('_rsl', '')] = gdfi['MLLW'] - gdfi[scen] # if this is 1, then we lost it (1-0)
        gdfi['tile'] = path.stem.rstrip('_beach')
        gdfo = gdfi[f"{scen.replace('_rsl', '')} geometry tile poly_ix".split()]
        lst_gsers.append(gdfo)
        
    gdf_reg = pd.concat(lst_gsers)
    gdf_reg.to_file(dst)
    log.critical ('Wrote:', dst)
    return


if __name__ == '__main__':
    region = 'South'
    habit  = 'rocky'
    scen   = 'med_rsl2050'
    path_wd = Path(os.getenv('dataroot')) / 'Sea_Level' / 'SFEI'
    log.critical (f'Begun {region}\n')
    main(region, habit, scen, path_wd, use_s2=False, test=False)
    # concat_beaches(region, scen, path_wd, use_s2=False)
