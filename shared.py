""" 
Global variables
"""
import os 
from pathlib import Path
import os.path as op
from dataclasses import dataclass
from typing import Literal
import warnings
import re

import time
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as xrr

from osgeo import gdal
from shapely.geometry import Point
from fiona.drvsupport import supported_drivers
supported_drivers['LIBKML'] = 'rw'
gdal.SetConfigOption('OGR_GEOJSON_MAX_OBJ_SIZE', '64000MB') # for fiona loading

warnings.filterwarnings("ignore", message="The nodata value *", category=UserWarning, module="rioxarray.raster_writer")


def df2gdf(path_df, attrcols=[], dst=None, epsg=4326):
    """
    Convert the df to a geo dataframe or projected (depending on epsg)
    Save as a point shpfile if dst is specified
    Driver can be "ESRI Shapefile" or GeoJSON

    Col names of lon/lat columns OR xt, default epgs is wgs84 (lat/lon)
    attrcol can be one string or list of strings
    """
    if isinstance(path_df, gpd.GeoDataFrame):
        print ('Already a geodataframe')
        return path_df

    elif isinstance(path_df, (pd.DataFrame, pd.Series)):
        df   = path_df
    else:
        df   = pd.read_pickle(path_df)
        dst  = op.join(op.dirname(path_df), op.splitext(op.basename(path_df))[0])

    df         = df.copy()

    if df.dropna(how='all').empty: 
        print('Dataframe is empty or full of Nans')
        return gpd.GeoDataFrame()

    if int(epsg) == 4326:
        raise Exception('Need to implement get_lalo_cols')
        latc, lonc, T, geoc = [*get_lalo_cols(df), 'lalo']
        if T: df = df.T

    elif int(epsg) in [4087, 3411, 3717, 26911]: # new PlateCarree
        latc, lonc, geoc = 'y', 'x', 'xy'
    else:
        raise bbError(f'EPSG:{epsg} is not currently supported')

    if isinstance(attrcols, str) and not str == 'all':
        attrcols = df.columns.tolist() if attrcols == 'all' else [attrcols]

    attrcols.remove(lonc) if lonc in attrcols else ''
    attrcols.remove(latc) if latc in attrcols else ''

    df[geoc] = df.apply(lambda df: Point((float(df[lonc]), float(df[latc]))), axis=1)

    cols   = [geoc] if not attrcols else [geoc, *attrcols]
    df     = df[cols]
    gdf    = gpd.GeoDataFrame(df, geometry=geoc, crs=f'EPSG:{epsg}')
    # print (geo_df.to_crs(epsg='4087') ## also identical
    if dst is not None:
        dst = f'{dst}.GeoJSON' if not dst.endswith('GeoJSON') else dst
        gdf.to_file(dst)
        print (f'Wrote geodataframe to:\n\t{dst}')
    return gdf


class CARIRegion(object):
    """
    base class for CARI region (North, Central, South)
    Sets up the paths to the working director yand the DEM
    """

    def __init__(self, region, path_wd=None, use_s2=True):
        assert region.title() in 'South Central North'.split(), 'Choose from South / Central / North'
        self.region = region.title()
        self.path_wd: str | Path = Path(os.getenv('dataroot')) / 'Sea_Level' / 'SFEI'  if path_wd is None else path_wd
        self.use_s2 = use_s2
        self.path_res = Path(self.path_wd) / 'results'
        self.setup_region()

    
    def setup_region(self):
        path_dems = self.path_wd / 'dems'
        dct_reg_info = {'Central': ['CA_Central_CoNED_DEM_2017_8657', 'Cent*[A-Za-z][0-9]', 3717, 'ca2017_central_coned_m8657.kmz'],
                        'North': ['CA_north_coned_DEM_2020_9181', f'Northern*_[0-9][0-9][0-9]', 3717, 
                                  'tileindex_CA_north_coned_DEM_2020.shp'],
                        'South': ['CA_Southern_CoNED_DEM_2016_8658', f'Southern*_[0-9][0-9]', 26911, 'socal.kmz']
                       }
        stem_ca_dems, self.stem, self.epsg, stem_bounds = dct_reg_info[self.region]
        
        self.path_ca_dems = path_dems / stem_ca_dems
        gdf_coned_bounds = gpd.read_file(self.path_ca_dems / stem_bounds).to_crs(self.epsg)
        self.Wr, self.Sr, self.Er, self.Nr = gdf_coned_bounds.total_bounds
        return

    
    def get_cari(self, kind='beaches'):
        assert kind.lower() in 'beaches rocky+beaches rocky'.split(), \
        'Incorrect CARI type, choose "beaches", "rocky+beaches", or "rocky"'
        if kind == 'beaches':
            path_cari = 'S2_Class_Polys/FINAL BEACH VECTOR.shp' if self.use_s2 else 'CARI_polygons/Cari_Beach.shp'
        elif kind == 'rocky':
            path_cari = 'CARI_polygons/Cari_Rocky.shp'
        else:
            raise Exception('Not implemented')
        
        gdf_cari = gpd.read_file(self.path_wd / path_cari)
        gdf_cari_coned = gdf_cari.to_crs(self.epsg)
        gdf_cari_coned.cx[self.Wr:self.Er, self.Sr:self.Nr]
        return gdf_cari_coned


class SetupProj(CARIRegion):
    def __init__(self, region, kind='beach', scen0='med_rsl2050', path_wd=None, use_s2=False):
        super().__init__(region, path_wd, use_s2)
        self.scen0 = scen0
        self.kind0 = kind.lower()
        self.set_slr('MLLW')
        self.set_slr('MHW')
        self.gdf_cari0 = self.get_cari(self.kind0)

    
    def set_slr(self, kind='MLLW'):
        """ These are made in MLLW_SLR.ipynb """
        kind = kind.upper()
        assert kind in 'MLLW MHW'.split(), 'Choose MLLW or MHW'
        lst_das = []
        for scen in f'0 {self.scen0}'.split():
            da_mllw = xrr.open_rasterio(self.path_wd / f'{kind}_SLR_{scen}.tif')
            da_mllw_re = da_mllw.sel(band=1).rio.reproject(self.epsg)
            da_mllw_re = da_mllw_re.where(da_mllw_re < 1e20, np.nan)
            da_mllw_re.rio.write_nodata(da_mllw_re.rio.nodata, encoded=True, inplace=True)

            # crop it to the coned region (N, Central, S)
            da_mllw_re_crop = da_mllw_re.sel(x=slice(self.Wr, self.Er), y=slice(self.Nr, self.Sr))
            lst_das.append(da_mllw_re_crop.assign_attrs(scenario=scen))

        if kind == 'MLLW':
            self.da_mllw0, self.da_mllw_slr = lst_das
        elif kind == 'MHW':
            self.da_mhw0, self.da_mhw_slr = lst_das
            
        return