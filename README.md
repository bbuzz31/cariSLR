# cariSLR

Hazard assessment of SLR to CA Coast using CARI: https://www.sfei.org/cari

## Install
Use conda/mamba/pip to install the python packages in `requirements.txt`

## Required Files
- CA State SLR Scenarios (included here)
- CARI Geodatabase: https://www.sfei.org/data
- Interpolated MLLW/MHHW+MAH: derived and converted into NAVD88 using NOAA's VDatum tool (https://vdatum.noaa.gov/). Contact SFEI for details.
- InSAR VLM (optional): https://zenodo.org/records/11154177

## Required Directory Structure
- The algorithms allow the specification of a working directory (`WD`)
- The CARI polygons must be located at `WD/`
- the MLLW/MHHW+MAH files must be located at `WD/JPL_Share/`
- (optional): CA_VLM.tif must be located at `WD/`

## Overview
- `RunTest.ipynb` contains a brief overview of the routines and a test for a single polygon
- `make_slr_scenario.py` will add the 5 SLR scenario (including VLM if specified) to interpolated and NAVD-referenced MLLW and MAH layers provided by SFEI
- `make_slr_cari.py` is the actual script that compares the DEMs to the SLR scenarios for each polygon and concatenates the results
- the additional notebooks are for understanding edgecases by viewing intermediate results as well as computing and visualizing final statistics 

 
