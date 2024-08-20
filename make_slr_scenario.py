from shared import *
import logging

def get_tidal_datum(path_wd, kind='MLLW'):
    kind = kind.upper()
    assert kind in 'MLLW MHW MAH'.split(), 'Choose MLLW, MHW, or MAH'
    kind = 'MHHW_Plus_MAH' if kind == 'MAH' else kind
    
    path_td = path_wd / 'JPL_Share' / f'{kind}_NAVD88m_IDW.tif'
    ds = xrr.open_rasterio(path_td)
    da = ds.sel(band=1)
    da = da.where(da < 1e20, np.nan)
    da.rio.write_nodata(da.rio.nodata, encoded=True, inplace=True)
    return da


def get_scenario_new(path_wd, scenario='Int', quant=50, year='2050'):
    df = pd.read_excel(path_wd.parent / 'Projections' / 'CA' / 'CA_SW__scenarios.xlsx')
    ser = df[((df['scenario'] == scenario) & (df['quantile'] == quant))]
    slr = ser[int(year)].item() # units = mm
    return slr/1e3 


def get_vlm(path_wd, overwrite=False):
    from shapely.geometry import box
    da0 = get_tidal_datum(path_wd, 'MLLW') # 300 m
    dst = path_wd / 'CA_VLM_MLLW.tif'
    if dst.exists() and not overwrite:
        print ('Using existing VLM')
        da = xrr.open_rasterio(dst).sel(band=1)
        da = da.where(da < 1e20)
        return da.rio.write_nodata(da.rio.nodata, encoded=True, inplace=True)
        
    src = path_wd / 'CA_VLM.tif'
    da_vlm = xrr.open_rasterio(src).sel(band=1).assign_attrs(units='mm/yr')
    da_vlm_re = da_vlm.rio.reproject(da0.rio.crs)
    da_vlm_re_clip = da_vlm_re.rio.clip([box(*da0.rio.bounds())], all_touched=True)
    da_vlm_re_in = da_vlm_re_clip.interp_like(da0, method='linear')
    da_vlm_re_in = da_vlm_re_in.where(da_vlm_re_in < 1e20)
    da_vlm_re_in.rio.write_nodata(da_vlm_re_in.rio.nodata, encoded=True, inplace=True)
    # doesnt really save nodata or units
    da_vlm_re_in.rio.to_raster(dst)
    print (f'Wrote reprojected, cropped, interpolated VLM to: {dst}')
    return da_vlm_re_in


def main(path_wd, year=2050, vlm=True):
    da_vlm = get_vlm(path_wd, overwrite=False)
    # convert to meters in the future and put 0 where nan
    da_vlm_proj = da_vlm * (year-2024) / 1000
    da_vlm_proj = da_vlm_proj.where(da_vlm_proj.notnull(), 0)
    for kind in 'MLLW MAH'.split():
        da0 = get_tidal_datum(path_wd, kind) # 300 m
        for scen in '0 Low IntLow Int IntHigh High'.split():
            fname = f'{kind}_SLR_{scen}{year}'
            if scen == '0':
                fname = fname.replace(str(year), '')
                slr = 0
                da_vlm = 0
            else:
                slr = get_scenario_new(path_wd, scenario=scen, year=year)
                da_vlm = da_vlm_proj
                
            da = da0 + slr 
            
            if vlm:
                da1 = da.copy()
                # fails on slr = 0
                try:
                    da1.data = da.data - da_vlm.data # x/y alignment
                except:
                    pass
                del da
                da = da1.copy()
                fname = f'{fname}_VLM'
                
            dst = path_wd / 'tidal_datums' / f'{fname}.tif'
            da.rio.to_raster(dst)
            log.info('Wrote: %s', dst)
            
if __name__ == '__main__':
    path_wd = Path('/scratch/tws_grace/data/Sea_Level/SFEI')
    main(path_wd, 2050)
    main(path_wd, 2100)
