import xarray as xr
import rioxarray as xrr
from shared import *

def compare_elevations_poly(poly, da_dem_in, da_mw_0, da_mw_slr, kind='MLLW'):
    """ The actual work: compare the DEM elevations within the poly against  current/future SLR """
    assert kind.upper() in 'MLLW MHW MTL'.split(), f'Incorrect kind: {kind}'
    if kind.upper() == 'MLLW':
        da_mw_0_poly = da_mw_0.rio.clip([poly], da_mw_0.rio.crs, all_touched=True)
        da_mw_slr_poly = da_mw_slr.rio.clip([poly], da_mw_slr.rio.crs, all_touched=True)
        
        da_dem_mw_0_poly = (da_dem_in>da_mw_0_poly).astype(int)
        da_dem_mw_slr_poly = (da_dem_in>da_mw_slr_poly).astype(int)
 
        # preserve nodata values in dem for more accurate percentages
        da_dem_mw_0_poly = da_dem_mw_0_poly.where(~da_dem_in.isnull())
        da_dem_mw_slr_poly = da_dem_mw_slr_poly.where(~da_dem_in.isnull())
        
        pct0  = 100*(da_dem_mw_0_poly.mean())
        pct1  = 100*(da_dem_mw_slr_poly.mean())
        pct_lost = pct0-pct1
        msg = f'\n\tLoss due to MLLW: {pct0.item():.3f} {pct1.item():.3f}. Lost: {pct_lost.item():.3f}% ' 
        print (msg)
        
        scen = da_mw_slr.attrs['scenario']
        df_mw = da_dem_mw_0_poly.rename('MLLW').to_dataframe()['MLLW'].reset_index().dropna()
        df_mw_slr = da_dem_mw_slr_poly.rename(scen).to_dataframe()[scen].reset_index().dropna()
        df_mw[scen] = df_mw_slr[scen]
    return df_mw


def make_rocky_beaches(region, scen='med_rsl2050', path_wd=None, use_s2=False, test=False):
    tstl     = '_test' if test else ''
    s2       = '_s2' if use_s2 else ''
    assert not use_s2, 'Sentinel2 not supported yet'
    kind0    = 'rocky'
    Obj      = SetupProj(region, kind0, scen, path_wd, use_s2)
    path_log = Obj.path_res / f'log_{Obj.region.title()}_{kind0}_pct{tstl}{s2}.txt'
    logf     = open(path_log, 'a')
    st0      = time.time()
    gdf_enso_map = get_enso_map(Obj.epsg) # enso dem path and its bounds
    
    lst_gsers = []
    # iterate over all DEM tiles
    for i, dem in enumerate((Obj.path_ca_dems).glob(f'{Obj.stem}.tif')):
        path_res = dem.parent / f'{dem.stem}_{kind0}{tstl}{s2}.GeoJSON'
        # dem = Obj.path_ca_dems / 'CentCA_south_Topobathy_CoNED_1m_B2.tif' for testing
        if path_res.exists() and not test:
            print (f'gdf of beaches in tile {dem.stem} exists, skipping...')
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
            msg = f'No {kind0} polygons within {dem.stem}'
            print (msg)
            print (msg, file=logf)
            continue
        else:
            print (f'{n_polys} {kind0} polygons within {dem.stem}')
        
        # check there is actual coned data within the polygons
        da_dem_check = da_dem.rio.clip(gdf_cari_tile.geometry.tolist(), gdf_cari_tile.crs, 
                                       drop=True, invert=False, all_touched=True)
        if da_dem_check.isnull().all():
            msg  = f'No actual DEM values within ANY {kind0} polygons. Skipping {dem.stem}.'
            print (msg)
            print (msg, file=logf)
            continue
            
        
        # interpolate the mllw and mhw to the 1m coned tiles
        da_mllw_re     = Obj.da_mllw0.interp_like(da_dem, method='linear')
        da_mllw_slr_re = Obj.da_mllw_slr.interp_like(da_dem, method='linear')
        
        da_mhw_re     = Obj.da_mhw0.interp_like(da_dem, method='linear')
        da_mhw_slr_re = Obj.da_mhw_slr.interp_like(da_dem, method='linear')
        
        elapi = time.time()-sti
        msg   = f'Preparing DEM, MLLW, MHW took {elapi/60:.2f} min.'
        print (msg)
        print (msg, file=logf)
        logf.flush()

        # iterate over each polygon in this tile
            ## this could probably be vectorized, as in da_dem_check
        lst_gdfs = []
        for j, (poly_ix, row) in enumerate(gdf_cari_tile.iterrows()):
            stj = time.time()
            poly = row.geometry
            cari_id = row['CARI_id']
            
            da_dem1 = da_dem.rio.clip([poly], gdf_cari_tile.crs, drop=True, invert=False, all_touched=True)
            if da_dem1.isnull().all():
                msg = f'No actual DEM values within {kind0} polygon {poly_ix} of {dem.stem}, skipping this polygon.'
                print (msg)
                print (msg, file=logf)
                continue
        
            # now get the mllw in this polygon
            da_mllw0i    = da_mllw_re.rio.clip([poly], gdf_cari_tile.crs, drop=True, invert=False, all_touched=True)
            da_mllw_slri = da_mllw_slr_re.rio.clip([poly], gdf_cari_tile.crs, drop=True, invert=False, all_touched=True)
            
            # get the mhw in this  polygon
            da_mhw0i    = da_mhw_re.rio.clip([poly], gdf_cari_tile.crs, drop=True, invert=False, all_touched=True)
            da_mhw_slri = da_mhw_slr_re.rio.clip([poly], gdf_cari_tile.crs, drop=True, invert=False, all_touched=True)

        
            # ------- compare elevations against MLLW (coned criteria 1) --- ELEVATION ABOVE MLLW
            msg = f'{dem.stem}, polyix={poly_ix} cari id={cari_id}: '\
            print (msg)
            print (msg, file=logf)
            df_mllw = compare_elevations_poly(poly, da_dem1, da_mllw0i, da_mllw_slri, 'MLLW')

            # ------ compare CONED elevations against MHW --- ELEVATION BELOW MHW
            da_dem_mhw0i    = (da_dem1<da_mhw0i).astype(int)
            da_dem_mhw_slri = (da_dem1<da_mhw_slri).astype(int)
            # preserve nodata values in dem for more accurate percentages
            da_dem_mhw0i    = da_dem_mhw0i.where(~da_dem1.isnull())
            da_dem_mhw_slri = da_dem_mhw_slri.where(~da_dem1.isnull())
            pct0_mhw  = 100*(da_dem_mhw0i.mean())
            pct1_mhw  = 100*(da_dem_mhw_slri.mean())
            pct_lost_mhw = pct0_mhw-pct1_mhw
            
            # ------- compare ENSO elevations against MHW 
            gdf_enso_map_poly = gdf_enso_map[gdf_enso_map.contains(poly)]
            if gdf_enso_map_poly.empty:
                msg= (f'No ENSO dem contains polygon cari_id={cari_id}, poly_ix={poly_ix}, skipping.') 
                print (msg)
                print (msg, file=logf) 
                continue
                
            elif gdf_enso_map_poly.shape[0]>1:
                msg = (f'More than 1 ENSO dem contains polygon cari_id={cari_id}, poly_ix={poly_ix}, skipping (for now).') 
                print (msg)
                print (msg, file=logf) 
                continue
                
            da_enso_dem = xrr.open_rasterio(gdf_enso_map_poly['path'].item())
            da_enso_dem_re = da_enso_dem.rio.reproject(Obj.epsg).sel(band=1)
            da_enso_dem_re = da_enso_dem_re.where(da_enso_dem_re>-1e20)
            da_enso_dem_re.rio.write_nodata(da_enso_dem_re.rio.nodata, encoded=True, inplace=True)
            # downsample from 0.5 m to 1 m horizontal resolution to match CoNED and MLLW/MHW
            da_enso_dem_re1 = da_enso_dem_re.interp_like(da_dem1, method='linear')
            da_enso_dem_re_poly = da_enso_dem_re1.rio.clip([poly], da_enso_dem_re1.rio.crs, drop=True, invert=False, all_touched=True)
            
            da_enso_mhw0i    = (da_enso_dem_re_poly<da_mhw0i).astype(int)
            da_enso_mhw_slri = (da_enso_dem_re_poly<da_mhw_slri).astype(int)
            
            # preserve nodata values in dem for more accurate percentages
            da_enso_mhw0i    = da_enso_mhw0i.where(~da_enso_dem_re_poly.isnull())
            da_enso_mhw_slri = da_enso_mhw_slri.where(~da_enso_dem_re_poly.isnull())
            pct0_mhw_enso  = 100*(da_enso_mhw0i.mean())
            pct1_mhw_enso  = 100*(da_enso_mhw_slri.mean())
            pct_lost_mhw_enso = pct0_mhw_enso-pct1_mhw_enso
            
            
            # msg = 
                  # f'\n\tLoss due to Coned MHW: {pct0_mhw.item():.3f} {pct1_mhw.item():.3f}. Lost: {pct_lost_mhw.item():.3f}% ' \
                  # f'\n\tLoss due to ENSO MHW: {pct0_mhw_enso.item():.3f} {pct1_mhw_enso.item():.3f}. Lost: {pct_lost_mhw_enso.item():.3f}%'
                  # f'\n\tLoss due to MLLW: {pct0.item():.3f} {pct1.item():.3f}. Lost: {pct_lost.item():.3f}% ' \
            
            # print (msg)
            # print (msg, file=logf)
            logf.flush()
            
            scen    = da_mllw_slr_re.attrs['scenario']
            
            df_mllw['poly_ix'] = poly_ix

            
            df_mhw = da_dem_mhw0i.rename('MHW').to_dataframe()['MHW'].reset_index().dropna()
            df_mhw_slr = da_dem_mhw_slri.rename(scen).to_dataframe()[scen].reset_index().dropna()
            df_mhw[f'{scen}_MHW'] = df_mhw_slr[scen]

            ## ENSO
            df_mhw_enso = da_enso_mhw0i.rename('MHW').to_dataframe()['MHW'].reset_index().dropna()
            df_mhw_slr_enso = da_enso_mhw_slri.rename(scen).to_dataframe()[scen].reset_index().dropna()
            df_mhw_enso[f'{scen}_MHW_ENSO'] = df_mhw_slr_enso[scen]

            ## at least for the coned, we should have all the same locations, so the same sizes.
            df_m = pd.merge(df_mllw, df_mhw, on='x y'.split())
            df_m = pd.merge(df_m, df_mhw_enso, on='x y'.split())
            assert df_m.shape[0] == df_mllw.shape[0] == df_mhw.shape[0] == df_mhw_enso.shape[0], 'Shapes arent the same?'
            
            lst_gdfs.append(df2gdf(df_m, 'all', epsg=Obj.epsg))

            if test and j == 0:
                return da_dem1, da_mllw0i, da_mllw_slri, df_m, poly
                
        gdf_tile = pd.concat(lst_gdfs)
        gdf_tile.to_file(path_res)
        
        msg = f'Wrote: {path_res}'
        print (msg)
        print (msg, file=logf)

    elap_tot = time.time()-st0
    msg = (f'Finished writing epsg={Obj.epsg} CoNED Rocky Percents in {elap_tot/60:.2f} min.')
    print (msg)
    print (msg, file=logf)
    print (f'Log file at:', path_log)
    return 


def concat_beaches(region, scen='med_rsl2050', path_wd=None, use_s2=False):
    """ Concatenate the beaches into a single (big) geojson file (e.g. 'South_beaches.GeoJSON') """
    Obj = CARIRegions(region, path_wd, scen, use_s2)
    s2_ext = '_s2' if use_s2 else ''
    paths_res = sorted(list(Obj.path_ca_dems.glob(f'*beach{s2_ext}.GeoJSON')))
    dst = Obj.path_res / f'{region}_beaches{s2_ext}.GeoJSON' 

    print (f'Concatenating: {len(paths_res)} {region} CoNED tiles')
    lst_gsers = []
    for i, path in enumerate(paths_res):
        print (f'Tile {i} of {len(paths_res)}') if i % 5 == 0  else ''
        gdfi  = gpd.read_file(path)
        gdfi[scen.replace('_rsl', '')] = gdfi['MLLW'] - gdfi[scen] # if this is 1, then we lost it (1-0)
        gdfi['tile'] = path.stem.rstrip('_beach')
        gdfo = gdfi[f"{scen.replace('_rsl', '')} geometry tile poly_ix".split()]
        lst_gsers.append(gdfo)
        
    gdf_reg = pd.concat(lst_gsers)
    gdf_reg.to_file(dst)
    print ('Wrote:', dst)
    return


if __name__ == '__main__':
    region = 'South'
    scen   = 'med_rsl2050'
    path_wd = Path(os.getenv('dataroot')) / 'Sea_Level' / 'SFEI'
    
    print (f'Begun {region}')
    make_rocky_beaches(region, scen, path_wd, use_s2=False, test=True)
    # concat_beaches(region, scen, path_wd, use_s2=False)
