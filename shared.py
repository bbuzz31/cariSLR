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
import xarray as xr
import rioxarray as xrr

from osgeo import gdal
from shapely.geometry import Point, box
from fiona.drvsupport import supported_drivers
try:
    from cariLog import log2file, logger as log
except:
    pass

supported_drivers['LIBKML'] = 'rw'
gdal.SetConfigOption('OGR_GEOJSON_MAX_OBJ_SIZE', '64000MB') # for fiona loading

warnings.filterwarnings("ignore", message="The nodata value *", 
                        category=UserWarning, module="rioxarray.raster_writer")
#### --------------------------------------------------------------------------------- Globals
regions = 'South Central North'.split()
# scen = scen0.replace('_rsl', '')
# s2_ext = '_s2' if use_s2 else ''
# vlm_ext = '_vlm' if use_vlm else ''

try:
    import contextily as cx
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.io import img_tiles as cimgt
    
    basemap_d = cimgt.GoogleTiles(url='https://server.arcgisonline.com/ArcGIS/rest/services/Elevation/World_Hillshade_Dark/MapServer/tile/{z}/{y}/{x}.jpg')
    basemap   = cimgt.GoogleTiles(url='https://server.arcgisonline.com/arcgis/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}.jpg')
    # basemap   = cimgt.GoogleTiles(url='https://www.google.cn/maps/vt?lyrs=s@189&gl=cn&x={x}&y={y}&z={z}'))
    
    cxbasemap   = cx.providers.Esri.WorldImagery
    cxbasemap_d = cx.providers.CartoDB.DarkMatter
    cxbasemap_t = cx.providers.USGS.USTopo
    projp = ccrs.PlateCarree()
    dct_proju = {'Central': ccrs.UTM(10), 'North': ccrs.UTM(10), 'South': ccrs.UTM(11)}
    
except:
    print ('Cannot plot. Probably cartopy/contextily not installed')
    
GFS = GS  = 14 
TFS = 20 # title
XFS = YFS = 16 # x/y axis
GFS = 14 # deg marks
CFS = 18 # colorbar labels

pt_parms  = dict(verticalalignment='top', fontsize=20, color='w',
                bbox=dict(boxstyle='square,pad=0.05', facecolor='black'))


sty_cari = dict(facecolor='none', edgecolor='deeppink') 
sty_coned = dict(facecolor='none', edgecolor='red', linestyle='--')


#### --------------------------------------------------------------------------------- Helper Functions

def df2gdf(path_df, attrcols=[], dst=None, epsg=4326):
    """
    Convert the df to a geo dataframe or projected (depending on epsg)
    Save as a point shpfile if dst is specified
    Driver can be "ESRI Shapefile" or GeoJSON

    Col names of lon/lat columns OR xt, default epgs is wgs84 (lat/lon)
    attrcol can be one string or list of strings
    """
    if isinstance(path_df, gpd.GeoDataFrame):
        log.debug ('Already a geodataframe')
        return path_df

    elif isinstance(path_df, (pd.DataFrame, pd.Series)):
        df   = path_df
    else:
        df   = pd.read_pickle(path_df)
        dst  = op.join(op.dirname(path_df), op.splitext(op.basename(path_df))[0])

    df         = df.copy()

    if df.dropna(how='all').empty: 
        log.error('Dataframe is empty or full of Nans')
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
        log.info (f'Wrote geodataframe to:\n\t{dst}')
    return gdf


def cartopy_cbar(im, xlabel='', orient='vertical', fontsize=13, labelsize=13, **plot_kws):
    # if passing an axes, try to get children; doesnt work yet
    # im     = im.get_children()[0] if 'GeoAxes' in type(im) else im
    im = im.get_children()[0] if isinstance(im, cartopy.mpl.geoaxes.GeoAxesSubplot) else im
    dct_ps = dict(ylabel='', size='5%', pad=0.15, ticks=None, pack_start=False)
    dct_ps.update(**plot_kws)

    divider = make_axes_locatable(im.axes)
    if orient == 'vertical':
        cbar_ax = divider.new_horizontal(size=dct_ps['size'], pad=dct_ps['pad'],
                            axes_class=plt.Axes, pack_start=dct_ps['pack_start'])
    else:
        cbar_ax = divider.new_vertical(size=dct_ps['size'], pad=dct_ps['pad'],
                            axes_class=plt.Axes)#, pack_start=dct_ps['pack_start'])

    cbar_ax = im.get_figure().add_axes(cbar_ax)
    cbar_ax.set_xlabel(xlabel) if isinstance(xlabel, str) else cbar_ax.set_xlabel(xlabel[0], size=xlabel[1], labelpad=10)
    cbar_ax.set_yscale('linear')
    cbar_ax.tick_params(labelsize=labelsize)

    cbar    = plt.colorbar(im, cax=cbar_ax, ticks=dct_ps['ticks'], spacing='uniform', extend=dct_ps.get('extend', 'neither'))
    cbar.set_label(dct_ps['ylabel'], rotation=270, labelpad=dct_ps.get('labelpad', 20), fontsize=fontsize)

    return cbar


def savefigs(path, figures=False, overwrite=False, extend='', alpha=0, **plot_kws):
    """
    Save current fig(s) to path
    Default save to path/Figures
    Optionally overwrite
    Optionally add an extension, e.g. in a loop for different grids
    Optionally set background transparent (default)
    plot_kws can contain dpi or fmt
    """
    path_figdir = op.join(path, 'Figures') if figures else path
    ext_fmt     = plot_kws.get('fmt', mpl.rcParams['savefig.format'])
    if not op.isdir(path_figdir):
        os.makedirs(path_figdir)
    figs = list(map(plt.figure, plt.get_fignums()))
    for fig in figs:

        fig.patch.set_alpha(alpha)

        #mpl.rcParams['figure.figsize']   = (18, 12) # figure out a way to make this work
        label    = fig.get_label() if fig.get_label() else 'tmp'
        path_fig = op.join(path_figdir, f'{label}.{ext_fmt}')
        i = 1
        while op.exists(path_fig):
            if overwrite: break
            path, ext = op.splitext(path_fig)
            ## overwrite last number
            if re.search(r'\d+$', path) is not None:
                path = path.rsplit('-', 1)[0]
            path_fig = f'{path}-{i}{ext}'
            label    = op.splitext(op.basename(path_fig))[0]
            i += 1

        if len(extend) > 0: label += str(extend)
        label_clean = label.replace('.', '-')
        dpi = plot_kws.get('dpi', fig.dpi)
        kws = dict(dpi=dpi, bbox_inches='tight', pad_inches=0.025, transparent=False,
                    edgecolor=fig.get_edgecolor())
        kws.update(plot_kws)

        dst = op.join(path_figdir, f'{label_clean}.{ext_fmt}')
        fig.savefig(dst, **kws)
        log.info('Saved figure: %s', dst)
        if 'SSH_CONNECTION' in os.environ:
            plt.close(fig)


def get_enso_map(path_enso, epsgi, verbose=False):
    """ make a map of ENSO filenames to their bounds in a specific EPSG """
    from shapely.geometry import box
    paths, polys = [], []
    dst = path_enso / f'enso_map_{epsgi}.GeoJSON'
    if dst.exists():
        # print (f'Using existing enso filename map for EPSG: {epsgi}') if verbose else ''
        return gpd.read_file(dst)

    lst_dems = sorted(list((path_enso / 'UTM10').glob('*.tif')))
    if epsgi == 26911 or epsgi == 4326:
        lst_dems = sorted(lst_dems + list((path_enso / 'UTM11').glob('*.tif')))
    
    ## ENSO is UTM10, want to convert to 3717 or 26911(South) (or 4326)
    for i, enso in enumerate(lst_dems):
        da_enso = xrr.open_rasterio(enso).sel(band=1)
        gser_enso = gpd.GeoSeries(box(*da_enso.rio.bounds()), crs=da_enso.rio.crs)
        poly = gser_enso.to_crs(epsgi).geometry.item()
        if i % 100 == 0:
            log.debug (f'Projecting bounds to: {enso.stem}, {i} of {len(lst_dems)}')

        paths.append(str(enso))
        polys.append(poly)

    gdf = gpd.GeoDataFrame(paths, geometry=polys, columns=['path'], crs=epsgi)
    gdf.to_file(dst)
    log.info (f'Wrote enso to dem map for EPSG: {epsgi}')
    return gdf
    

#### --------------------------------------------------------------------------------- Classes

class CARIRegion(object):
    """
    base class for CARI region (North, Central, South)
    Sets up the paths to the working director yand the DEM
    """

    def __init__(self, region, path_wd=None, use_s2=False):
        assert region.title() in 'South Central North'.split(), 'Choose from South / Central / North'
        self.region = region.title()
        self.path_wd: str | Path = Path(os.getenv('dataroot')) / 'Sea_Level' / 'SFEI'  if path_wd is None else path_wd
        self.use_s2 = use_s2
        self.path_res = Path(self.path_wd) / 'results'
        self.path_enso_dems = self.path_res.parent / 'dems' / 'West_Coast_El_Nino'
        self.path_res.mkdir(exist_ok=True)
        self.path_enso_dems.mkdir(exist_ok=True)
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
        layer = 'kml_image_ca2017_central_coned_m8657_L1_0_0' if self.region == 'Central' else None
        gdf_coned_bounds = gpd.read_file(self.path_ca_dems / stem_bounds, layer=layer).to_crs(self.epsg)
        self.Wr, self.Sr, self.Er, self.Nr = gdf_coned_bounds.total_bounds
        return

    
    def get_cari(self, kind='beach', new=True):
        assert kind.lower() in 'beach rocky rocky_mllw'.split(), \
        'Incorrect CARI type, choose "beach", or "rocky"'
        ext = 'GeoJSON' if new else 'shp'
        log.debug(f'Using new polygons? {new}')
        if kind == 'beach':
            path_cari = 'S2_Class_Polys/FINAL BEACH VECTOR.shp' if self.use_s2 else f'CARI_polygons/Cari_Beach.{ext}'
        elif kind in 'rocky rocky_mllw'.split():
            path_cari = 'S2_Class_Polys/FINAL ROCKY VECTOR.shp' if self.use_s2 else f'CARI_polygons/Cari_Rocky.{ext}'
        else:
            raise Exception('Not implemented')
        
        gdf_cari = gpd.read_file(self.path_wd / path_cari)
        gdf_cari_coned = gdf_cari.to_crs(self.epsg)
        gdf_cari_coned.cx[self.Wr:self.Er, self.Sr:self.Nr]
        if self.use_s2:
            gdf_cari_coned.drop(columns='id path'.split(), inplace=True)
            gdf_cari_coned['CARI_id'] = gdf_cari_coned.index # dummy for later on
        return gdf_cari_coned


class SetupProj(CARIRegion):
    def __init__(self, region, habit='beach', scen0='Int2050', path_wd=None, use_vlm=False, use_s2=False):
        super().__init__(region, path_wd, use_s2)
        self.scen0 = scen0
        self.habit0 = habit.lower()
        self.use_vlm = use_vlm
        self.set_slr('MLLW')
        self.set_slr('MAH')
        self.gdf_cari0 = self.get_cari(self.habit0)

    
    def set_slr(self, kind='MLLW'):
        """ These are made in MLLW_SLR.ipynb """
        vlm = '_VLM' if self.use_vlm else ''
        kind = kind.upper()
        assert kind in 'MLLW MHW MAH'.split(), 'Choose MLLW, MHW, or MAH'
        lst_das = []
        for scen in f'0 {self.scen0}'.split():
            da_mllw = xrr.open_rasterio(self.path_wd / 'tidal_datums' / f'{kind}_SLR_{scen}{vlm}.tif')
            da_mllw_re = da_mllw.sel(band=1).rio.reproject(self.epsg)
            da_mllw_re = da_mllw_re.where(da_mllw_re < 1e20, np.nan)
            da_mllw_re.rio.write_nodata(da_mllw_re.rio.nodata, encoded=True, inplace=True)

            # crop it to the coned region (N, Central, S)
            da_mllw_re_crop = da_mllw_re.sel(x=slice(self.Wr, self.Er), y=slice(self.Nr, self.Sr))
            lst_das.append(da_mllw_re_crop.assign_attrs(scenario=scen))

        if kind == 'MLLW':
            self.da_mllw_0, self.da_mllw_slr = lst_das
        elif kind == 'MHW':
            self.da_mhw_0, self.da_mhw_slr = lst_das
        elif kind == 'MAH':
            self.da_mah_0, self.da_mah_slr = lst_das
 
        return
