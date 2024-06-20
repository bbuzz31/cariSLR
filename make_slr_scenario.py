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


def main(path_wd, year=2050):
    for kind in 'MLLW MAH'.split():
        da0 = get_tidal_datum(path_wd, kind)
        for scen in '0 Low IntLow Int IntHigh High'.split():
            fname = f'{kind}_SLR_{scen}{year}'
            if scen == '0':
                fname = fname.replace(str(year), '')
                slr = 0
            else:
                slr = get_scenario_new(path_wd, scenario=scen, year=year)
            da = da0 + slr
            dst = path_wd / 'tidal_datums' / f'{fname}.tif'
            da.rio.to_raster(dst)
            log.info('Wrote: %s', dst)
            
if __name__ == '__main__':
    path_wd = Path('/scratch/tws_grace/data/Sea_Level/SFEI')
    main(path_wd, 2050)
    main(path_wd, 2100)
