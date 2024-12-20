import xarray as xr
from rioxarray.merge import merge_arrays
import logging
import datetime
from shared import *


def get_enso_poly(gdf_enso_map_in, da_dem_in, poly_in):
    """ Get the ENSO dem values within a polygon and reproject/interpolate to the CoNED """
    gdf_enso_map_poly = gdf_enso_map_in[gdf_enso_map_in.intersects(poly_in)]
    if gdf_enso_map_poly.empty:
        return 0

    lst_da_enso = []
    for path in gdf_enso_map_poly['path']:
        da_enso_dem1 = xrr.open_rasterio(path)
        lst_da_enso.append(da_enso_dem1)
        
    da_enso_dem = merge_arrays(lst_da_enso) if len(lst_da_enso) > 1 else lst_da_enso[0]
    
    da_enso_dem_re = da_enso_dem.rio.reproject(da_dem_in.rio.crs).sel(band=1)
    da_enso_dem_re = da_enso_dem_re.where(da_enso_dem_re>-1e20)
    da_enso_dem_re.rio.write_nodata(da_enso_dem_re.rio.nodata, encoded=True, inplace=True)
    # downsample from 0.5 m to 1 m horizontal resolution to match CoNED and MLLW/MAH
    da_enso_dem_re1 = da_enso_dem_re.interp_like(da_dem_in, method='linear')
    
    # try:
    # except Exception as E:
    #     log.error(E)
    #     da_enso_poly = path
    # lst_da_enso_polys.append(da_enso_poly)
    # lst_da_enso_polys.append(da_enso_dem_re1)

    
    # clip to polygon
    # ## this has to do because it SOUTH and needs other
    
    return da_enso_dem_re1.rio.clip([poly_in], da_dem_in.rio.crs, all_touched=True)
            

def compare_elevations_poly(poly_in, da_dem_in, da_mw_0_in, da_mw_slr_in, kind='MLLW'):
    """ The actual work: compare the DEM elevations within the poly against current/future SLR 

    For MLLW (beach), each pixel has a 1 if it's above mllw (dry, good) or below mllw (always wet, bad)
    For MAH (rocky), each pixel has a 1 if it's below MAH (good, rocky) or 0 if it's above (always dry)
    """
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
        da_dem_mw_0_poly = da_dem_mw_0_poly.where(da_dem_in.notnull())
        da_dem_mw_slr_poly = da_dem_mw_slr_poly.where(da_dem_in.notnull())
        
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


def main(region, habit='rocky', scen='Int2050', path_wd=None, use_vlm=False, use_s2=False, test=False):
    """ Run the actual analysis; Save estimates of losses per each DEM in the dems folder as GeoJSON """
    habit = habit.lower()
    tstl = '_test' if test else ''
    s2 = '_s2' if use_s2 else ''
    vlm = '_vlm' if use_vlm else ''
    today = datetime.datetime.today().date().strftime('%Y%m%d')
    assert habit in 'beach rocky rocky_mllw'.split(), f'Incorrect habit: {habit}'

    Obj  = SetupProj(region, habit, scen, path_wd, use_vlm, use_s2)
    path_log = Obj.path_res / f'log_{Obj.region.title()}_{habit}_pct{tstl}{vlm}{s2}_{today}.txt'
    log2file(path_log)
    st0  = time.time()
    gdf_enso_map = get_enso_map(Obj.path_enso_dems, Obj.epsg) # enso dem path and its bounds
    
    lst_gsers = []
    # iterate over all DEM tiles
    for i, dem in enumerate((Obj.path_ca_dems).glob(f'{Obj.stem}.tif')):
        path_res = dem.parent / f'{dem.stem}_{habit}{tstl}{vlm}{s2}_{scen}.GeoJSON'
        # dem = Obj.path_ca_dems / 'CentCA_south_Topobathy_CoNED_1m_B3.tif' # for testing
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
            print (''); log.info(f'{n_polys} {habit} polygons within {dem.stem}')
        
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
            
            da_dem1 = da_dem.rio.clip([poly], gdf_cari_tile.crs, 
                                      drop=True, invert=False, all_touched=True).assign_attrs(tile=dem.stem)
            try:
                da_dem.rio.clip([poly], gdf_cari_tile.crs, all_touched=False)
            except:
                log.warning(f'The DEM only lightly touches polygon: {poly_ix} of {dem.stem}, skipping.')
                continue
                
            log.critical (f'{dem.stem}, polyix={poly_ix} cari id={cari_id}:')
            df_mllw = compare_elevations_poly(poly, da_dem1, da_mllw_re, da_mllw_slr_re, 'MLLW')
            
            if habit == 'rocky':
                df_mah = compare_elevations_poly(poly, da_dem1, da_mah_re, da_mah_slr_re, 'MAH')
                da_enso_dem = get_enso_poly(gdf_enso_map, da_dem1, poly)
                if isinstance(da_enso_dem, int):
                    log.warning (f'{da_enso_dem} ENSO dems with polygon cari_id={cari_id}, poly_ix={poly_ix}, skipping.') 
                    continue
                    
                df_mah_enso = compare_elevations_poly(poly, da_enso_dem, da_mah_re, da_mah_slr_re, 'MAH_ENSO')
                df_m = pd.merge(df_mllw, df_mah, on='x y'.split(), how='outer')
                df_m = pd.merge(df_m, df_mah_enso, on='x y'.split(), how='outer')
                del df_mah, da_enso_dem, df_mah_enso
            else:
                df_m = df_mllw
            df_m['poly_ix'] = poly_ix
            df_m['cari_id'] = cari_id
            
            ## enso may be smaller...
            # assert df_m.shape[0] == df_mllw.shape[0] == df_mah.shape[0] == df_mah_enso.shape[0], 'Shapes arent the same?'
            
            lst_gdfs.append(df2gdf(df_m, 'all', epsg=Obj.epsg))

            if test and (df_m['MLLW'] - df_m[f'{scen}_MLLW']).sum() > 0:
                return da_dem1, da_mllw_re, da_mllw_slr_re, df_m, poly

            del df_m, da_dem1, df_mllw

        if not lst_gdfs:
            log.warning(f'Polygons in CONED tile {dem.stem} had no ENSO values, skipping.')
            continue
            
        gdf_tile = pd.concat(lst_gdfs)
        gdf_tile.to_file(path_res)
        log.critical(f'CONED tile completed: Wrote: {path_res}')
        del gdf_tile, lst_gdfs, da_dem, da_dem_check, gdf_cari_tile

    elap_tot = time.time()-st0
    log.critical (f'Finished writing {region} {habit} Percents in {elap_tot/60:.2f} min.')
    log.info (f'Log file at: %s', path_log)
    return 


def concat_results(region, habit, scen='Int2050', path_wd=None, use_vlm=False, use_s2=False):
    """ Concatenate the beaches/rocky geojson results into a single csv """
    habit = habit.lower()
    assert habit in 'beach rocky rocky_mllw'.split(), f'Incorrect habit: {habit}'
    Obj = SetupProj(region, habit, scen, path_wd, use_vlm, use_s2)
    s2_ext = '_s2' if use_s2 else ''
    vlm_ext = '_vlm' if use_vlm else ''
    paths_res = sorted(list(Obj.path_ca_dems.glob(f'*{habit}{vlm_ext}{s2_ext}_{scen}.GeoJSON')))
    
    if not paths_res:
        log.warning(f'No results found for {habit}{vlm_ext}{s2_ext}_{scen}')
        return
    
    dst = Obj.path_res / f'{region}_{habit}_{scen}{vlm_ext}{s2_ext}.csv' 

    log.info (f'Concatenating: {len(paths_res)} {region} {habit} CoNED tiles')
    lst_gsers = []
    col0 = f'{scen}_MLLW'
    for i, path in enumerate(paths_res):
        log.info (f'Tile {i} of {len(paths_res)}') if i % 5 == 0  else ''
        gdfi  = gpd.read_file(path)
        # if this is 1, then we lost it; 1 is above MLLW 0 is below (1-0)
        # - overwrites column name
        gdfi[f'{scen}_MLLW'] = gdfi['MLLW'] - gdfi[f'{scen}_MLLW'] 
        gdfi['tile'] = path.stem.rstrip(f'_{habit}')
        cols = col0
        
        if habit == 'rocky':
            # temporary, for sanity 
            gdfi = gdfi.rename(columns={'MAH': 'MAH_CONED', f'{scen}_MAH': f'{scen}_MAH_CONED'})
            # if either is 1, then we call it good and rocky
            # note that because of the OR we may still gain on one condition but since it was already rocky on other condition there will be no change
            gdfi['MAH'] = np.where((gdfi['MAH_CONED'] + gdfi['MAH_ENSO']) > 0, 1, 0) 
            gdfi[f'{scen}_MAH'] = np.where((gdfi[f'{scen}_MAH_CONED'] + gdfi[f'{scen}_MAH_ENSO']) > 0, 1, 0) 
            # if these 1, we lost it, if 0 no change, if -1, we gained
            gdfi[f'{scen}_MAH'] = gdfi['MAH'] - gdfi[f'{scen}_MAH']  # -- this is final result

            # these overwrite and are kept in case want to look individually
            gdfi[f"{scen}_MAH_CONED"] = gdfi['MAH_CONED'] - gdfi[f'{scen}_MAH_CONED']
            gdfi[f"{scen}_MAH_ENSO"] = gdfi['MAH_ENSO'] - gdfi[f'{scen}_MAH_ENSO']
            cols = f'{col0} {scen}_MAH {scen}_MAH_CONED {scen}_MAH_ENSO'
            
        gdfo = gdfi[f'{cols} poly_ix tile geometry'.split()]
        lst_gsers.append(gdfo)

    gdf_reg = pd.concat(lst_gsers)
    df_reg = pd.DataFrame(gdf_reg) # smaller and easier
    df_reg.to_csv(dst)
    log.critical ('Wrote: %s', dst)
    return


def concat_results_poly(habit='beach', years=[2050, 2100], path_wd=None, use_vlm=False, use_s2=False):
    """ Calculate the area/percent lost in each polygon and write to GeoJSON """
    path_wd = Path(os.getenv('dataroot')) / 'Sea_Level' / 'SFEI'  if path_wd is None else path_wd
    path_res = path_wd / 'results'
    s2_ext = '_s2' if use_s2 else ''
    vlm_ext = '_vlm' if use_vlm else ''
    scens = 'Low IntLow Int IntHigh High'.split()
    years = [years] if not isinstance(years, list) else years
    # scens = 'High'.split()
    log.info(f'Calculating average within {habit} polygons')
    for i, yeari in enumerate(years):
        lst_res = []
        for sceni in scens:
            sy = f'{sceni}{yeari}'
            lst_paths = sorted(list(path_res.glob(f'*_{habit}_{sy}{vlm_ext}{s2_ext}.csv')))
    
            for path in lst_paths:
                reg, hab, sce =  path.stem.split('_')[:3]
                dfi = pd.read_csv(path, index_col=0)
                # do this in here cuz of different crs
                BObj = CARIRegion(reg, path_wd, use_s2)
                gdf_cari = BObj.get_cari(hab).to_crs(4326)
                for poly, dfj in dfi.groupby('poly_ix'):
                    tile = dfj['tile'].unique()
                    # 1 means we lost it (MLLW) or (MAH)
                    # -1 means we gained it (MAH)
                    pct_lost = 100 * dfj[f'{sy}_MLLW'].mean()
                    geom = gdf_cari.loc[poly].geometry
                    if hab == 'beach':
                        lst_res.append((poly, pct_lost, sceni, yeari, reg, geom)) 
                        cols = 'cari_ix pct_lost scenario year region geometry'.split()
                    elif hab == 'rocky':
                        # these are below MLLW in future; they are technically 'gained' cuz 
                        # z was greater than MAH (too high) but now it's below; but they are below MLLW which matters more
                        df_wrong = dfj[(dfj[f'{sy}_MLLW'] + dfj[f'{sy}_MAH'].abs()) > 1]
                        dfj.loc[df_wrong.index, f'{sy}_MAH'] = 0
                        pct_gain = -100 * dfj[f'{sy}_MAH'].mean()
                        
                        col_total = (dfj[f'{sy}_MAH'] + dfj[f'{sy}_MLLW'])
                        # this doesnt happen
                        if col_total.max() > 1:
                            print ('How are we losing from both MLLW and MAH?')
                            breakpoint()
                        pct_total = -100*col_total.mean()
                        lst_res.append((poly, pct_lost, pct_gain, pct_total, sceni, yeari, reg, geom)) 
        cols = 'cari_ix pct_lost pct_gain pct_total scenario year region geometry'.split()
        df_res_poly = pd.DataFrame(lst_res, columns=cols)
        gdf_res_poly = gpd.GeoDataFrame(df_res_poly)
        
        ## get rid of the duplicate polygons (partly in different regions) and replace with mean
        dupe_polys = []
        new_rows = []
        for ix, dfi in gdf_res_poly.groupby('scenario year cari_ix'.split()):
            if dfi.shape[0] > 1:
                row = dfi.iloc[0].copy()
                row['pct_lost'] = dfi['pct_lost'].mean()
                if habit == 'rocky':
                    row['pct_gain'] = dfi['pct_gain'].mean()
                    row['pct_total'] = dfi['pct_total'].mean()
                dupe_polys.append(row['cari_ix'].item())
                new_rows.append(row)
        
        gdf_res_de = gdf_res_poly[~gdf_res_poly['cari_ix'].isin(dupe_polys)]
        if dupe_polys:
            # gdf_res_de = pd.concat([gdf_res_de, pd.concat(new_rows)]).sort_index()
            gdf_res_de = pd.concat([gdf_res_de, gpd.GeoDataFrame(new_rows)]).sort_index()

        # NAD83 / California Albers, units = m 
        gdf_res_de['area'] = gdf_res_de.set_crs(4326).to_crs(3310).area
        gdf_res_de['area_lost'] = gdf_res_de['area'] * (gdf_res_de['pct_lost']*0.01)
        if hab == 'rocky':
            gdf_res_de['area_gain'] = gdf_res_de['area'] * (gdf_res_de['pct_gain']*0.01)

        dst = path_res / f'{hab}_polygons_lost_{yeari}{vlm_ext}{s2_ext}.GeoJSON' 
        gdf_res_de.to_file(dst)
        log.info(f'Wrote: {dst}')
    return


if __name__ == '__main__':
    habits  = 'rocky'.split()
    use_vlm = False
    use_s2  = False
    path_wd = Path(os.getenv('dataroot')) / 'Sea_Level' / 'SFEI'
    for habit in habits:
        # for scen0 in 'Low IntLow Int IntHigh High'.split():
        for scen0 in 'High'.split():
            for year in [2050]:#, 2100]:
                for region in 'Central'.split():
                    scen = f'{scen0}{year}'
                    # log.critical (f'Begun {habit} {region}, {scen}\n')
                    main(region, habit, scen, path_wd, use_s2=use_s2, use_vlm=use_vlm, test=False)
                    concat_results(region, habit, scen, path_wd, use_vlm=use_vlm, use_s2=use_s2)


    ## only when all scenarios/years are done for CARI/CARI_VLM/S2/S2_VLM
    # concat_results_poly('beach', [2050, 2100], path_wd, use_vlm, use_s2)
    # concat_results_poly('rocky', [2050, 2100], path_wd, use_vlm, use_s2)

    # use_vlm=True
    # concat_results_poly('beach', [2050, 2100], path_wd, use_vlm, use_s2)
    # concat_results_poly('rocky', [2050, 2100], path_wd, use_vlm, use_s2)

