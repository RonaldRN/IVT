# 💧 Integrated Water Vapor Transport (IVT) over South America

This repository contains **Python scripts for calculating, analyzing, and visualizing Integrated Water Vapor Transport (IVT)** over South America, with particular emphasis on the Amazon Basin and its surrounding regions.

The project combines atmospheric fields from the **ERA5 reanalysis** with precipitation estimates from **IMERG** to investigate the climatological behavior and variability of atmospheric moisture transport.

The analyses include monthly and seasonal climatologies, precipitation patterns, ENSO-related anomalies, geopotential height anomalies, and long-term IVT trends.

---

# 🎯 Objectives

The main objectives of this repository are to:

* calculate vertically integrated atmospheric moisture transport from ERA5 data;
* characterize the monthly and seasonal climatology of IVT;
* analyze the zonal and meridional components of moisture transport;
* compare IVT patterns with precipitation from IMERG;
* investigate IVT and precipitation anomalies during selected **El Niño and La Niña events**;
* analyze geopotential height anomalies associated with these events;
* visualize spatial trends in IVT over the 2001–2020 period;
* and provide Python routines for atmospheric moisture-transport analysis over South America.

---

# 💧 Integrated Water Vapor Transport

The Integrated Water Vapor Transport combines information about **atmospheric moisture and horizontal wind** throughout a vertical atmospheric layer.

In vector form, the zonal and meridional components can be represented as:

[
IVT_u = -\frac{1}{g}\int q,u,dp
]

[
IVT_v = -\frac{1}{g}\int q,v,dp
]

where:

* (q) is the specific humidity;
* (u) is the zonal wind component;
* (v) is the meridional wind component;
* (p) is atmospheric pressure;
* (g) is gravitational acceleration.

The IVT magnitude is then calculated as:

[
IVT = \sqrt{IVT_u^2 + IVT_v^2}
]

In the current implementation of this repository, the vertical integration is performed between:

```text
1000–500 hPa
```

using the trapezoidal integration method.

The calculated products are stored as:

```text
IVT
IVT_u
IVT_v
```

allowing both the magnitude and direction of atmospheric moisture transport to be analyzed.

---

# 🛰️ Datasets

## ERA5

Atmospheric variables used for calculating IVT are obtained from the **ERA5 reanalysis**.

The main variables are:

```text
Specific humidity              q
Zonal wind component           u
Meridional wind component      v
```

ERA5 monthly-mean pressure-level data are used for the climatological and anomaly analyses.

The repository also includes routines for downloading:

* specific humidity;
* zonal wind;
* meridional wind;
* geopotential.

The IVT climatologies are primarily analyzed over the period:

```text
2001–2020
```

---

## 🌧️ IMERG

Precipitation analyses use the:

**IMERG — Integrated Multi-satellitE Retrievals for GPM**

product.

IMERG precipitation is used to construct:

* monthly climatologies;
* seasonal climatologies;
* precipitation anomalies;
* comparisons between precipitation and atmospheric moisture transport.

The main climatological period used in the current scripts is:

```text
2001–2020
```

---

# 🌎 Study region

The analyses focus on a large sector of South America, with particular emphasis on the **Amazon Basin and surrounding regions**.

Several plotting routines use approximately:

```text
Latitude:   40°S – 15°N
Longitude:  85°W – 30°W
```

allowing the moisture transport entering and leaving the Amazon region to be analyzed within its broader atmospheric context.

---

# 📂 Repository structure

The repository is currently organized as:

```text
IVT/
│
├── codigos/
│   │
│   ├── IMERG.py
│   ├── IVT.py
│   ├── IVT_Functions.py
│   ├── IVT_trend_2001-2020.py
│   ├── era5_download_data.py
│   ├── era5_download_data_altg.py
│   └── hgt.py
│
├── figuras/
│
└── README.md
```

---

# 🐍 Description of the scripts

## `IVT_Functions.py`

Contains the main functions used to calculate the vertically integrated moisture transport.

The module includes routines for:

* IVT calculations for specific periods;
* monthly IVT climatology;
* seasonal IVT climatology;
* zonal moisture transport (`IVT_u`);
* meridional moisture transport (`IVT_v`);
* total IVT magnitude.

This module is imported by other scripts in the repository.

---

## `IVT.py`

Main script for the **calculation and visualization of IVT**.

It uses the functions defined in:

```text
IVT_Functions.py
```

to produce analyses including:

* monthly IVT climatology;
* seasonal IVT climatology;
* IVT vectors;
* selected ENSO events;
* seasonal IVT anomalies.

Examples currently analyzed include:

### El Niño

```text
2015–2016
```

### La Niña

```text
2010–2011
```

The resulting maps combine IVT magnitude with the corresponding zonal and meridional transport components.

---

## `IMERG.py`

Contains routines for analyzing precipitation from IMERG.

The script produces:

* monthly precipitation climatology;
* seasonal precipitation climatology;
* precipitation anomalies for selected periods;
* El Niño precipitation anomalies;
* La Niña precipitation anomalies.

The precipitation climatology currently uses:

```text
2001–2020
```

This allows changes in moisture transport to be compared with changes in the spatial distribution of precipitation.

---

## `hgt.py`

Analyzes geopotential fields from ERA5 and produces **geopotential anomalies** associated with selected ENSO events.

The script constructs seasonal climatologies and compares individual events against the climatological reference.

Examples include:

```text
El Niño 2015–2016
La Niña 2010–2011
```

These fields provide large-scale atmospheric circulation context for the IVT and precipitation anomalies.

---

## `IVT_trend_2001-2020.py`

Produces maps of IVT trends for the period:

```text
2001–2020
```

using trend results obtained from the **Mann–Kendall statistical test**.

> **Important:** in the current version of the repository, this script reads a NetCDF file containing previously calculated Mann–Kendall trend results and produces the corresponding maps. The calculation of the Mann–Kendall test itself is not included in this script.

---

## `era5_download_data.py`

Uses the **Copernicus Climate Data Store (CDS) API** to obtain monthly ERA5 pressure-level data.

The requested atmospheric variables include:

```text
specific_humidity
u_component_of_wind
v_component_of_wind
```

These variables provide the atmospheric information required to calculate IVT.

---

## `era5_download_data_altg.py`

Uses the CDS API to download ERA5:

```text
geopotential
```

at:

```text
300 hPa
```

for the period:

```text
2001–2020
```

These data are used to investigate the large-scale circulation associated with the analyzed events.

---

# 🌊 ENSO analyses

The repository also investigates changes in atmospheric moisture transport and precipitation during selected phases of the **El Niño–Southern Oscillation (ENSO)**.

The workflow compares particular events against the 2001–2020 climatological reference.

Examples include:

```text
El Niño 2015–2016
        │
        ├── IVT anomaly
        ├── precipitation anomaly
        └── geopotential anomaly
```

and:

```text
La Niña 2010–2011
        │
        ├── IVT anomaly
        ├── precipitation anomaly
        └── geopotential anomaly
```

Additional figures for more recent ENSO periods are also currently available in the repository.

---

# 📈 Trend analysis

The repository includes IVT trend maps for:

```text
2001–2020
```

based on the non-parametric **Mann–Kendall trend test**.

The analysis is designed to identify regions where vertically integrated moisture transport exhibits systematic changes during the study period.

The plotting script uses previously calculated trend fields stored in:

```text
ivt_trend_2001-2020.nc
```

---

# 🖼️ Figures

The directory:

```text
figuras/
```

contains examples of products generated from the analyses.

These include:

```text
IVT monthly climatology
IVT seasonal climatology
IVT anomalies
IMERG monthly climatology
IMERG seasonal climatology
IMERG anomalies
geopotential anomalies
ENSO composites
IVT trends
```

Examples currently available include:

```text
ivt_clim.png
ivt_season.png
ivt_anom_El-Nino.png
ivt_anom_La-Nina.png
ivt_trend.png

imerg_monthly_clim.png
imerg_season_climatology.png
imerg_anom_El-Nino_2015-2016.png
imerg_anom_La-Nina_2010-2011.png

hgt_anom_El-Nino_2015-2016.png
hgt_anom_La-Nina_2010-2011.png
```

---

# 🔬 General workflow

The general computational workflow can be summarized as:

```text
ERA5 pressure-level data
(q, u, v)
        │
        ▼
Vertical integration
1000–500 hPa
        │
        ▼
┌─────────────────────────────┐
│            IVT              │
│                             │
│   IVT_u     IVT_v     IVT   │
└─────────────────────────────┘
        │
        ▼
Monthly and seasonal climatology
2001–2020
        │
        ├───────────────────────────┐
        │                           │
        ▼                           ▼
 ENSO anomalies             Long-term trends
        │                           │
        ▼                           ▼
El Niño / La Niña           Mann–Kendall
        │
        ▼
Comparison with IMERG
precipitation
        │
        ▼
Geopotential anomalies
        │
        ▼
Physical interpretation of
moisture transport variability
```

---

# ⚙️ Requirements

The scripts are written in **Python** and use packages from the scientific Python ecosystem.

Main dependencies include:

```text
numpy
pandas
xarray
matplotlib
cartopy
cmocean
cdsapi
netCDF4
```

A possible installation using Conda is:

```bash
conda install -c conda-forge numpy pandas xarray matplotlib cartopy cmocean netcdf4 cdsapi
```

or using `pip`:

```bash
pip install numpy pandas xarray matplotlib cartopy cmocean netCDF4 cdsapi
```

---

# 📥 Installation

Clone the repository:

```bash
git clone https://github.com/RonaldRN/IVT.git
```

Enter the repository:

```bash
cd IVT
```

and then:

```bash
cd codigos
```

---

# 🌐 ERA5 data access

To use the ERA5 download scripts, access to the **Copernicus Climate Data Store API** is required.

After configuring your CDS API credentials, ERA5 data can be downloaded using:

```bash
python era5_download_data.py
```

and:

```bash
python era5_download_data_altg.py
```

The requested years, variables, pressure levels and geographical domain can be modified directly in the corresponding scripts.

---

# ▶️ Running the analyses

The Python scripts can be executed from the command line.

For example:

```bash
python IVT.py
```

or:

```bash
python IMERG.py
```

Before running the analyses, verify the input paths and filenames defined inside each script.

---

# ⚠️ Reproducibility note

The current repository primarily contains the **analysis and plotting scripts** and selected output figures.

Several external input files referenced by the scripts are not currently distributed with the repository, including ERA5 and IMERG NetCDF datasets, topography data, climatological products and some geographic auxiliary files.

Examples include:

```text
AS_era5-montly_2001-2020.nc
ivt_climatology_monthly.nc
ivt_climatology_season.nc
imerg_monthly_2001-2020.nc
imerg_climatology_monthly_2001-2020.nc
imerg_climatology_season_2001-2020.nc
hgt_climatology_season_2001-2020.nc
ivt_trend_2001-2020.nc
topo_25.1.nc
```

Therefore, users interested in reproducing the complete workflow should obtain the required datasets and adapt the paths and filenames in the scripts to their local computational environment.

---

# 📌 Scientific interpretation

IVT provides a useful framework for investigating how atmospheric moisture is transported across South America.

Combining IVT with precipitation, geopotential fields and climate variability makes it possible to investigate questions such as:

* How does moisture transport vary throughout the annual cycle?
* Which regions act as important pathways for atmospheric moisture?
* How does moisture transport toward and across the Amazon Basin change seasonally?
* How are precipitation anomalies related to changes in moisture transport?
* How does atmospheric circulation differ between El Niño and La Niña conditions?
* Are there systematic changes in IVT over the 2001–2020 period?

---

# 👨‍💻 Author

**Ronald Guiuseppi Ramírez Nina**
Atmospheric Sciences
Institute of Astronomy, Geophysics and Atmospheric Sciences
University of São Paulo — IAG/USP

* 📧 [ronald.ramirez.nina@usp.br](mailto:ronald.ramirez.nina@usp.br)
* 📧 [ronald.ramirez.nina@alumni.usp.br](mailto:ronald.ramirez.nina@alumni.usp.br)
* 📧 [ronald.ramirez.nina@gmail.com](mailto:ronald.ramirez.nina@gmail.com)
* 🐙 [GitHub — @RonaldRN](https://github.com/RonaldRN)

---

# 📚 Citation

If you use or adapt the scripts available in this repository for scientific or academic purposes, please acknowledge the repository:

```text
RonaldRN/IVT
```

Repository:

https://github.com/RonaldRN/IVT

If these analyses are associated with a scientific publication, thesis, dissertation, or other research product, the corresponding citation can also be added here.

---

# 📄 License

A specific software license is not currently documented in this README.

If the repository is intended for reuse and redistribution by other researchers, adding a `LICENSE` file is recommended.

---

## 🌎 Acknowledgements

The analyses use atmospheric reanalysis data from **ERA5** and satellite precipitation estimates from **IMERG**.

The scripts were developed for atmospheric-science applications involving moisture transport, precipitation variability and large-scale circulation over South America.
 
