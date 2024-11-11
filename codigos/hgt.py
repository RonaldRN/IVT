'''
Código que cálcula as anomalias de altura geopotencial para períodos específicos
Faz o plot do ciclo anual e sazonal médio da precipitação
Faz o plot das anomalias de precipitação sazonal para Anos El Niño e La Niña
'''

# Importar livrarias
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.crs as ccrs  # Projeções do Cartopy
import cartopy.feature as cfeature  # Recursos geográficos
import cmocean
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.colors as mcolors
from IVT_Functions import IVT, IVT_clim_season, IVT_clim_monthly

# Carregar o arquivo .nc da hgt mensal
ds = xr.open_dataset('globe_era5_2001-2020_altg.nc')

# Converter a coordenada 'date' para datetime64
ds = ds.assign_coords(date=pd.to_datetime(ds['date'].astype(str), format='%Y%m%d'))
# Converter a longitude no formato de -180° até +180°
ds.coords['longitude'] = (ds.coords['longitude'] + 180) % 360 - 180
ds = ds.sortby(ds.longitude)

# Calcular a climatologia de DJF [2001-2020]   
#ds_hgt_season = ds.groupby('date.season').mean(dim=['date'])
# Salvar dado da média sazonal hgt
#ds_hgt_season.to_netcdf("hgt_climatology_season_2001-2020.nc") 

# Carregar o arquivo .nc da climatologia sazonal
ds_clim_hgt = xr.open_dataset('hgt_climatology_season_2001-2020.nc')
# Converter a longitude no formato de -180° até +180°
ds_clim_hgt.coords['longitude'] = (ds_clim_hgt.coords['longitude'] + 180) % 360 - 180
ds_clim_hgt = ds_clim_hgt.sortby(ds_clim_hgt.longitude)

# -------------------------------------------------------------------------
# Função para calcular anomalias
# -------------------------------------------------------------------------

def anomaly(ds, ds_clim, date1:str, date2:str, date3:str, season:str):

    # Selecionar a vari�vel 'q' (q-Specific humidity [kg. kg-1])
    ds_date1 = ds['z'].sel(date=date1)
    ds_date2 = ds['z'].sel(date=date2)
    ds_date3 = ds['z'].sel(date=date3)
    ds_mean = (ds_date1 + ds_date2 + ds_date3) / 3

    anomaly = ds_mean - ds_clim['z'].sel(season=season) 
    anomaly = anomaly/10

    return anomaly

# -------------------------------------------------------------------------
# Anomalias El Niño
# -------------------------------------------------------------------------

anomaly_nino_DJF = anomaly(ds, ds_clim_hgt, date1='2015-12-01', date2='2016-01-01', date3='2016-02-01', season='DJF')
anomaly_nino_MAM = anomaly(ds, ds_clim_hgt, date1='2016-03-01', date2='2016-04-01', date3='2016-05-01', season='MAM')
anomaly_nino_JJA = anomaly(ds, ds_clim_hgt, date1='2015-06-01', date2='2015-07-01', date3='2015-08-01', season='JJA')
anomaly_nino_SON = anomaly(ds, ds_clim_hgt, date1='2015-09-01', date2='2015-10-01', date3='2015-11-01', season='SON')

hgt_dataset_nino = [anomaly_nino_DJF, anomaly_nino_MAM, anomaly_nino_JJA, anomaly_nino_SON]

titles = ['Anomalia DJF - El Niño[2015-2016]\nhgt 300hPa', 'Anomalia MAM - El Niño[2016]\nhgt 300hPa',\
           'Anomalia JJA - El Niño[2015]\nhgt 300hPa', 'Anomalia SON - El Niño[2015]\nhgt 300hPa']

# Criar a figura com 1 linha e 2 colunas
fig, axes = plt.subplots(2, 2, figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})

# Ajustar para que axes seja uma lista plana para iterar facilmente
axes = axes.flatten()

# Loop para configurar e plotar cada subplot
for i, ax in enumerate(axes):
    # Adicionar linhas de costa, limites de países e rios
    ax.coastlines(resolution='110m', linewidth=1)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    
    # Agregando Estados
    states_provinces = cfeature.NaturalEarthFeature(
        category = 'cultural',
        name = 'admin_1_states_provinces_lines',
        scale = '50m',
        facecolor = 'none')
    ax.add_feature(states_provinces, edgecolor = '0.25')
    # Configurar as linhas de grade
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 12}
    gl.ylabel_style = {'size': 12}
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER

    # Configurar a norma para centralizar o branco entre -10 e 10
    norm = mcolors.TwoSlopeNorm(vmin=-120, vcenter=0, vmax=120)
    # Exemplo de plotagem do campo de magnitude do IVT para cada mês (ajustar conforme necessário)
    hgt_dataset_nino[i].plot(ax=ax, cmap=cmocean.cm.balance_r, norm=norm, transform=ccrs.PlateCarree(),
                             levels=np.arange(-120,121,30), vmin=-120, vmax=120, 
                             cbar_kwargs={'label':'[dmgp]', 'shrink':0.75, 'ticks':np.arange(-120, 121, 30)})
    
    levels = np.arange(-120, 121, 30)
    # Plotando contornos para valores negativos com linhas tracejadas
    neg_contours = ax.contour(ds_clim_hgt.longitude, ds_clim_hgt.latitude, hgt_dataset_nino[i].isel(pressure_level=0),
                              colors='black', linestyles = '--', linewidths = 0.7, levels=levels[levels < 0], 
                              transform=ccrs.PlateCarree())
    # Plotando contornos para valores positivos com linhas contínuas
    pos_contours = ax.contour(ds_clim_hgt.longitude, ds_clim_hgt.latitude, hgt_dataset_nino[i].isel(pressure_level=0), 
                              colors = 'black', linestyles = '-', linewidths = 0.7, levels=levels[levels > 0], 
                              transform=ccrs.PlateCarree())

    # Título para cada subplot (mês)
    ax.set_title(titles[i], fontsize=14)

# Ajustar o layout para evitar sobreposição
plt.subplots_adjust(hspace=0.15, wspace=0.15) #, left=0.125, bottom=0.1, right=0.9, top=0.9)
plt.tight_layout()
plt.show()
fig.savefig('hgt_anom_El-Nino_2015-2016.png', format='png', dpi = 300)

# -------------------------------------------------------------------------
# Anomalias El Niño
# -------------------------------------------------------------------------

anomaly_nina_DJF = anomaly(ds, ds_clim_hgt, date1='2010-12-01', date2='2011-01-01', date3='2011-02-01', season='DJF')
anomaly_nina_MAM = anomaly(ds, ds_clim_hgt, date1='2011-03-01', date2='2011-04-01', date3='2011-05-01', season='MAM')
anomaly_nina_JJA = anomaly(ds, ds_clim_hgt, date1='2010-06-01', date2='2010-07-01', date3='2010-08-01', season='JJA')
anomaly_nina_SON = anomaly(ds, ds_clim_hgt, date1='2010-09-01', date2='2010-10-01', date3='2010-11-01', season='SON')

hgt_dataset_nina = [anomaly_nina_DJF, anomaly_nina_MAM, anomaly_nina_JJA, anomaly_nina_SON]

titles = ['Anomalia DJF - La Niña[2010-2011]\nhgt 300hPa', 'Anomalia MAM - La Niña[2011]\nhgt 300hPa',\
           'Anomalia JJA - La Niña[2010]\nhgt 300hPa', 'Anomalia SON - La Niña[2010]\nhgt 300hPa']

# Criar a figura com 1 linha e 2 colunas
fig, axes = plt.subplots(2, 2, figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})

# Ajustar para que axes seja uma lista plana para iterar facilmente
axes = axes.flatten()

# Loop para configurar e plotar cada subplot
for i, ax in enumerate(axes):
    # Adicionar linhas de costa, limites de países e rios
    ax.coastlines(resolution='110m', linewidth=1)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    
    # Agregando Estados
    states_provinces = cfeature.NaturalEarthFeature(
        category = 'cultural',
        name = 'admin_1_states_provinces_lines',
        scale = '50m',
        facecolor = 'none')
    ax.add_feature(states_provinces, edgecolor = '0.25')
    # Configurar as linhas de grade
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 12}
    gl.ylabel_style = {'size': 12}
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER

    # Configurar a norma para centralizar o branco entre -10 e 10
    norm = mcolors.TwoSlopeNorm(vmin=-120, vcenter=0, vmax=120)
    # Exemplo de plotagem do campo de magnitude do IVT para cada mês (ajustar conforme necessário)
    hgt_dataset_nina[i].plot(ax=ax, cmap=cmocean.cm.balance_r, norm=norm, transform=ccrs.PlateCarree(),
                             levels=np.arange(-120,121,30), vmin=-120, vmax=120, 
                             cbar_kwargs={'label':'[dmgp]', 'shrink':0.75, 'ticks': np.arange(-120, 121, 30)})
    
    levels = np.arange(-120, 121, 30)
    # Plotando contornos para valores negativos com linhas tracejadas
    neg_contours = ax.contour(ds_clim_hgt.longitude, ds_clim_hgt.latitude, hgt_dataset_nina[i].isel(pressure_level=0),
                              colors='black', linestyles = '--', linewidths = 0.7, levels=levels[levels < 0], 
                              transform=ccrs.PlateCarree())
    # Plotando contornos para valores positivos com linhas contínuas
    pos_contours = ax.contour(ds_clim_hgt.longitude, ds_clim_hgt.latitude, hgt_dataset_nina[i].isel(pressure_level=0), 
                              colors = 'black', linestyles = '-', linewidths = 0.7, levels=levels[levels > 0], 
                              transform=ccrs.PlateCarree())

    # Título para cada subplot (mês)
    ax.set_title(titles[i], fontsize=14)

# Ajustar o layout para evitar sobreposição
plt.subplots_adjust(hspace=0.15, wspace=0.15) #, left=0.125, bottom=0.1, right=0.9, top=0.9)
plt.tight_layout()
plt.show()
fig.savefig('hgt_anom_La-Nina_2010-2011.png', format='png', dpi = 300)
