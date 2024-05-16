import xarray as xr
import rioxarray as xrr
from shared import *


class MakeBeaches(CARIRegion):
    def __init__(self, region, scen0='med_rsl2050', path_wd=None, use_s2=False):
        super().__init__(region, path_wd, use_s2)
        self.scen0 = scen0
        self.set_slr()
        self.gdf_beaches0 = self.get_beaches()

    
    def set_slr(self):
        """ These are made in MLLW_SLR.ipynb """
        lst_das = []
        for scen in f'0 {self.scen0}'.split():
            da_mllw = xrr.open_rasterio(self.path_wd / f'MLLW_SLR_{scen}.tif')
            da_mllw_re = da_mllw.sel(band=1).rio.reproject(self.epsg)
            da_mllw_re = da_mllw_re.where(da_mllw_re < 1e20, np.nan)
            da_mllw_re.rio.write_nodata(da_mllw_re.rio.nodata, encoded=True, inplace=True)

            # crop it to the coned region (N, Central, S)
            da_mllw_re_crop = da_mllw_re.sel(x=slice(self.Wr, self.Er), y=slice(self.Nr, self.Sr))
            lst_das.append(da_mllw_re_crop.assign_attrs(scenario=scen))
                                          
        self.da_mllw0, self.da_mllw_slr = lst_das
        return


def make_beaches(region, scen='med_rsl2050', path_wd=None, use_s2=False, test=False):
    tstl     = '_test' if test else ''
    s2       = '_s2' if use_s2 else ''
    Obj      = MakeBeaches(region, scen, path_wd, use_s2)
    path_log = Obj.path_wd / f'log_{Obj.region.title()}_pct{tstl}{s2}.txt'
    logf     = open(path_log, 'a')
    st0      = time.time()
    
    lst_gsers = []
    for i, dem in enumerate((Obj.path_ca_dems).glob(f'{Obj.stem}.tif')):
        path_res = dem.parent / f'{dem.stem}_beach{tstl}{s2}.GeoJSON'
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

        # crop beaches within this tile
        gdf_beach_tile = Obj.gdf_beaches0.cx[wi:ei, si:ni]
        n_polys = gdf_beach_tile.shape[0]
        if n_polys == 0:
            msg = f'No beach polygons within {dem.stem}'
            print (msg)
            print (msg, file=logf)
            continue
        else:
            print (f'{n_polys} beach polygons within {dem.stem}')
        
        # check there is actual coned data within the polygons
        da_dem_check = da_dem.rio.clip(gdf_beach_tile.geometry.tolist(), gdf_beach_tile.crs, 
                                       drop=True, invert=False, all_touched=True)
        if da_dem_check.isnull().all():
            msg  = f'No actual DEM values within ANY beach polygons. Skipping {dem.stem}.'
            print (msg)
            print (msg, file=logf)
            continue
            
        
        # interpolate the mllw to the 1m tiles
        da_mllw_re     = Obj.da_mllw0.interp_like(da_dem, method='linear')
        da_mllw_slr_re = Obj.da_mllw_slr.interp_like(da_dem, method='linear')
        
        elapi = time.time()-sti
        msg   = f'Preparing DEM and MLLW took {elapi/60:.2f} min.'
        print (msg)
        print (msg, file=logf)
        logf.flush()

        # iterate over each polygon in this tile
            ## this could probably be vectorized, as in da_dem_check
        lst_gdfs = []
        for j, (poly_ix, row) in enumerate(gdf_beach_tile.iterrows()):
        # for j, poly in enumerate(gdf_beach_tile.geometry):
            # poly_ix = gdf_beach_tile.index[j]
            stj = time.time()
            poly = row.geometry
            
            da_dem1 = da_dem.rio.clip([poly], gdf_beach_tile.crs, drop=True, invert=False, all_touched=True)
            if da_dem1.isnull().all():
                msg = f'No actual DEM values within beach polygon {poly_ix} of {dem.stem}, skipping this polygon.'
                print (msg)
                print (msg, file=logf)
                continue
        
            # now get the mllw in this polygon
            da_mllw0i    = da_mllw_re.rio.clip([poly], gdf_beach_tile.crs, drop=True, invert=False, all_touched=True)
            da_mllw_slri = da_mllw_slr_re.rio.clip([poly], gdf_beach_tile.crs, drop=True, invert=False, all_touched=True)
        
            # compare elevations
            da_dem_mllw0i    = (da_dem1>da_mllw0i).astype(int)
            da_dem_mllw_slri = (da_dem1>da_mllw_slri).astype(int)
            # preserve nodata values in dem for more accurate percentages
            da_dem_mllw0i    = da_dem_mllw0i.where(~da_dem1.isnull())
            da_dem_mllw_slri = da_dem_mllw_slri.where(~da_dem1.isnull())
            
            
            pct0  = 100*(da_dem_mllw0i.mean())
            pct1  = 100*(da_dem_mllw_slri.mean())
            pct_lost = pct0-pct1

            msg =  (f'{dem.stem}, poly ix{poly_ix}: {pct0.item():.3f} {pct1.item():.3f}. Lost: {pct_lost.item():.3f}%')
            print (msg)
            print (msg, file=logf)
            logf.flush()

            # do something with the results
            # if not pct0 == pct:
            scen    = da_mllw_slr_re.attrs['scenario']
            df_mllw = da_dem_mllw0i.rename('MLLW').to_dataframe()['MLLW'].reset_index().dropna()
            
            if df_mllw.empty:
                msg = f'No points left after dropping nans in {dem.stem} poly_ix {poly_ix}'
                print (msg)
                print (msg, file=logf)
                continue
                
            df_mllw_slr = da_dem_mllw_slri.rename(scen).to_dataframe()[scen].reset_index().dropna()
            df_mllw[scen] = df_mllw_slr[scen]
            df_mllw['poly_ix'] = poly_ix
            
            lst_gdfs.append(df2gdf(df_mllw, 'all', epsg=Obj.epsg))

            if test and j == 0:
                return da_dem1, da_mllw0i, da_mllw_slri, df_mllw, poly
                
        gdf_tile = pd.concat(lst_gdfs)
        gdf_tile.to_file(path_res)
        
        msg = f'Wrote: {path_res}'
        print (msg)
        print (msg, file=logf)

    elap_tot = time.time()-st0
    msg = (f'Finished writing epsg={Obj.epsg} CoNED Beach Percents in {elap_tot/60:.2f} min.')
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
    make_beaches(region, scen, path_wd, use_s2=False, test=False)
    concat_beaches(region, scen, path_wd, use_s2=False)
