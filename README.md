# cariSLR

Hazard assessment of SLR to CA Coast using CARI: https://www.sfei.org/cari

## Install
Use conda/mamba/pip to install the python packages in requirements.txt

## Overview

- `make_cari_polys.py` will create Cari_Beaches.GeoJSON and Cari_Rocky.GeoJSON from the geodatabase available from CARI
- `make_slr_scenario.py` will add the 5 SLR scenario (including VLM from Govorcin et al) to interpolated and NAVD-referenced MLLW and MAH layers provided by SFEI
- `make_slr_cari.py` is the actual script that compares the DEMs to the SLR scenarios for each polygon and concatenates the results
- the various notebooks are for understanding edgecases by viewing intermediate results as well as final statistics 
