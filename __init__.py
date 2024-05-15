""" 
Global variables
"""
import os 
from pathlib import Path
import os.path as op
from dataclasses import dataclass
from typing import Literal

import time
import numpy as np
import pandas as pd
import geopandas as gpd

from osgeo import gdal
from shapely.geometry import Point
from fiona.drvsupport import supported_drivers
supported_drivers['LIBKML'] = 'rw'
gdal.SetConfigOption('OGR_GEOJSON_MAX_OBJ_SIZE', '64000MB') # for fiona loading


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


class CARIRegion():
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

    
    def get_beaches(self):
        path_beach = 'S2_Class_Polys/FINAL BEACH VECTOR.shp' if self.use_s2 else 'CARI_polygons/Cari_Beach.shp'
        
        gdf_beach = gpd.read_file(self.path_wd / path_beach)
        gdf_beach_coned = gdf_beach.to_crs(self.epsg)
        gdf_beach_coned.cx[self.Wr:self.Er, self.Sr:self.Nr]
        return gdf_beach_coned
