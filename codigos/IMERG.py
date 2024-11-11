'''
Código que cálcula as anomalias de precipitação para períodos específicos
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

# Carregar o arquivo .nc da climatologia mensal
ds_clim_monthly = xr.open_dataset('imerg_climatology_monthly_2001-2020.nc')
ds_clim_monthly = ds_clim_monthly.sel(lat=slice(-40,15), lon=slice(-85,-30))

# Carregar o arquivo .nc da climatologia sazonal
ds_clim_season = xr.open_dataset('imerg_climatology_season_2001-2020.nc')
ds_clim_season = ds_clim_season.sel(lat=slice(-40,15), lon=slice(-85,-30))

# Carregar o arquivo .nc com os acumulados mensais do 2001-2020
ds = xr.open_dataset('imerg_monthly_2001-2020.nc')
ds = ds.sel(lat=slice(-40,15), lon=slice(-85,-30))

# Dataset da topografia
ds_topo = xr.open_dataset('topo_25.1.nc')
# Fazendo um slice para a �rea da Bacia Amazônica
ds_topo_ba = ds_topo.sel(lat=slice(-40, 10), lon=slice(-85, -40))

# ------------------------------------------------------------------------------------------------
# Figura da Climatologia [2001-2020] Mensal
# ------------------------------------------------------------------------------------------------
titles = ['Prec. Janeiro [2001-2020]', 'Prec. Fevereiro [2001-2020]', 'Prec. Março [2001-2020]',
          'Prec. Abril [2001-2020]', 'Prec. Maio [2001-2020]', 'Prec. Junho [2001-2020]', 'Prec. Julho [2001-2020]',
          'Prec. Agosto [2001-2020]', 'Prec. Setembro [2001-2020]', 'Prec. Outubro [2001-2020]',
          'Prec. Novembro [2001-2020]', 'Prec. Dezembro [2001-2020]']

# Número de dias dos mêses desde Janeiro até Dezembro
amount = [31,28,31,30,31,30,
          31,31,30,31,30,31]

# Criar a figura com 3 linha e 4 colunas
fig, axes = plt.subplots(3, 4, figsize=(15, 12), subplot_kw={'projection': ccrs.PlateCarree()})

# Ajustar para que axes seja uma lista plana para iterar facilmente
axes = axes.flatten()

# Loop para configurar e plotar cada subplot
for i, ax in enumerate(axes):
    # Adicionar linhas de costa, limites de países e rios
    ax.coastlines(resolution='110m', linewidth=1)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    
    # Exemplo de adição de rios
    rivers_50m = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '50m')
    ax.add_feature(rivers_50m, facecolor='None', edgecolor='b', linewidth=1)
    
    # Contorno da Bacia Amazônica
    fname = 'amazonas/Limite_Cuenca_Amazonas_geogpsperu_juansuyo.shp'
    adm1_shapes = list(shpreader.Reader(fname).geometries())
    ax.add_geometries(adm1_shapes, ccrs.PlateCarree(), edgecolor='yellow', facecolor='none', 
                      alpha=1, linewidth=1)
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

    # Exemplo de plotagem do campo de magnitude do IVT para cada mês (ajustar conforme necessário)
    # Unidades da precipitação são mm/hora por isso é necessário multiplicar por 24 e pelo número de dias do mês
    # para converter a mm/mês
    (ds_clim_monthly['prec'][i]*24*amount[i]).plot(ax=ax, cmap=cmocean.cm.rain, transform=ccrs.PlateCarree(),
                                   levels=np.arange(0,451,50), vmax=450, cbar_kwargs={'label':'[$mm.mes^{-1}$]', 'shrink':0.5})

    # Plotar a topografia
    colorb = np.arange(1500, 1550, 50)
    ax.contour(ds_topo_ba.lon, ds_topo_ba.lat, ds_topo_ba['z'], 
               levels=colorb, colors=['brown'], linewidths=0.5)

    # Título para cada subplot (mês)
    ax.set_title(titles[i], fontsize=14)

# Ajustar o layout para evitar sobreposi��o
plt.subplots_adjust(hspace=0.5, wspace=0.15) #, left=0.125, bottom=0.1, right=0.9, top=0.9)
plt.tight_layout()
plt.show()
fig.savefig('imerg_monthly_clim.png', format='png', dpi = 300)

# -------------------------------------------------------------------------------------------------
# Figura da Climatologia [2001-2020] Sazonal 
# ------------------------------------------------------------------------------------------------

titles = ['Prec. DJF [2001-2020]', 'Prec. MAM [2001-2020]', 'Prec. JJA [2001-2020]', 'Prec. SON [2001-2020]' ]

# Criar a figura com 2 linha e 2 colunas
fig, axes = plt.subplots(2, 2, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

# Ajustar para que axes seja uma lista plana para iterar facilmente
axes = axes.flatten()

season = ['DJF','MAM','JJA','SON']
# Loop para configurar e plotar cada subplot
for i, ax in enumerate(axes):
    # Adicionar linhas de costa, limites de países e rios
    ax.coastlines(resolution='110m', linewidth=1)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    
    # Exemplo de adição de rios
    rivers_50m = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '50m')
    ax.add_feature(rivers_50m, facecolor='None', edgecolor='b', linewidth=1)
    
    # Contorno da Bacia Amazônica
    fname = 'amazonas/Limite_Cuenca_Amazonas_geogpsperu_juansuyo.shp'
    adm1_shapes = list(shpreader.Reader(fname).geometries())
    ax.add_geometries(adm1_shapes, ccrs.PlateCarree(), edgecolor='yellow', facecolor='none', 
                      alpha=1, linewidth=1)
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

    # Exemplo de plotagem do campo de magnitude do IVT para cada mês (ajustar conforme necessário)
    (ds_clim_season['prec']*24*30).sel(season=season[i]).plot(ax=ax, cmap=cmocean.cm.rain, transform=ccrs.PlateCarree(),
                                   levels=np.arange(0,451,50), vmax=450, cbar_kwargs={'label':'[$mm.mes^{-1}$]', 'shrink':0.5})

    # Plotar a topografia
    colorb = np.arange(1500, 1550, 50)
    ax.contour(ds_topo_ba.lon, ds_topo_ba.lat, ds_topo_ba['z'], 
               levels=colorb, colors=['brown'], linewidths=0.5)

    # Título para cada subplot (mês)
    ax.set_title(titles[i], fontsize=14)

# Ajustar o layout para evitar sobreposi��o
plt.subplots_adjust(hspace=0.45, wspace=0.15) #, left=0.125, bottom=0.1, right=0.9, top=0.9)
plt.tight_layout()
plt.show()
fig.savefig('imerg_season_climatology.png', format='png', dpi = 300)

# -------------------------------------------------------------------------
# Função para calcular anomalias
# -------------------------------------------------------------------------

def anomaly(ds, ds_clim, date1:str, date2:str, date3:str, season:str):

    # Selecionar a vari�vel 'q' (q-Specific humidity [kg. kg-1])
    ds_date1 = ds['prec'].sel(time=date1)
    ds_date2 = ds['prec'].sel(time=date2)
    ds_date3 = ds['prec'].sel(time=date3)
    ds_mean = (ds_date1 + ds_date2 + ds_date3) / 3

    anomaly = ds_mean - ds_clim['prec'].sel(season=season) 
    anomaly = anomaly * 24 * 30         # Para converter de mm/hora -> mm/mês

    return anomaly

# -------------------------------------------------------------------------
# Anomalias El Niño
# -------------------------------------------------------------------------

anomaly_nino_DJF = anomaly(ds, ds_clim_season, date1='2015-12-01', date2='2016-01-01', date3='2016-02-01', season='DJF')
anomaly_nino_MAM = anomaly(ds, ds_clim_season, date1='2016-03-01', date2='2016-04-01', date3='2016-05-01', season='MAM')
anomaly_nino_JJA = anomaly(ds, ds_clim_season, date1='2015-06-01', date2='2015-07-01', date3='2015-08-01', season='JJA')
anomaly_nino_SON = anomaly(ds, ds_clim_season, date1='2015-09-01', date2='2015-10-01', date3='2015-11-01', season='SON')

anomaly_dataset_nino = [anomaly_nino_DJF, anomaly_nino_MAM, anomaly_nino_JJA, anomaly_nino_SON]

titles = ['Prec. Anomalia DJF - El Niño[2015-2016]', 'Prec. Anomalia MAM - El Niño[2016]',\
           'Prec. Anomalia JJA - El Niño[2015]', 'Prec. Anomalia SON - El Niño[2015]']

# Criar a figura com 1 linha e 2 colunas
fig, axes = plt.subplots(2, 2, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

# Ajustar para que axes seja uma lista plana para iterar facilmente
axes = axes.flatten()

# Loop para configurar e plotar cada subplot
for i, ax in enumerate(axes):
    # Adicionar linhas de costa, limites de países e rios
    ax.coastlines(resolution='110m', linewidth=1)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    
    # Exemplo de adição de rios
    rivers_50m = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '50m')
    ax.add_feature(rivers_50m, facecolor='None', edgecolor='b', linewidth=1)
    
    # Contorno da Bacia Amazônica
    fname = 'amazonas/Limite_Cuenca_Amazonas_geogpsperu_juansuyo.shp'
    adm1_shapes = list(shpreader.Reader(fname).geometries())
    ax.add_geometries(adm1_shapes, ccrs.PlateCarree(), edgecolor='yellow', facecolor='none', 
                      alpha=1, linewidth=1)
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
    norm = mcolors.TwoSlopeNorm(vmin=-100, vcenter=0, vmax=100)
    # Exemplo de plotagem do campo de magnitude do IVT para cada mês (ajustar conforme necessário)
    anomaly_dataset_nino[i].plot(ax=ax, cmap=cmocean.cm.tarn, norm=norm, transform=ccrs.PlateCarree(),
                            levels=np.arange(-100,101,10), vmin=-100, vmax=100, 
                            cbar_kwargs={'label':'[$mm.mes^{-1}$]', 'ticks':np.arange(-100, 101, 20), 'shrink':0.75})

    # Plotar a topografia
    colorb = np.arange(1500, 1550, 50)
    ax.contour(ds_topo_ba.lon, ds_topo_ba.lat, ds_topo_ba['z'], 
               levels=colorb, colors=['brown'], linewidths=0.5)

    # Título para cada subplot (mês)
    ax.set_title(titles[i], fontsize=14)

# Ajustar o layout para evitar sobreposi��o
plt.subplots_adjust(hspace=0.45, wspace=0.15) #, left=0.125, bottom=0.1, right=0.9, top=0.9)
plt.tight_layout()
plt.show()
fig.savefig('imerg_anom_El-Nino_2015-2016.png', format='png', dpi = 300)

# -------------------------------------------------------------------------
# Anomalias La Niña
# -------------------------------------------------------------------------

anomaly_nina_DJF = anomaly(ds, ds_clim_season, date1='2010-12-01', date2='2011-01-01', date3='2011-02-01', season='DJF')
anomaly_nina_MAM = anomaly(ds, ds_clim_season, date1='2011-03-01', date2='2011-04-01', date3='2011-05-01', season='MAM')
anomaly_nina_JJA = anomaly(ds, ds_clim_season, date1='2010-06-01', date2='2010-07-01', date3='2010-08-01', season='JJA')
anomaly_nina_SON = anomaly(ds, ds_clim_season, date1='2010-09-01', date2='2010-10-01', date3='2010-11-01', season='SON')

anomaly_dataset_nina = [anomaly_nina_DJF, anomaly_nina_MAM, anomaly_nina_JJA, anomaly_nina_SON]

titles = ['Prec. Anomalia DJF - La Niña[2010-2011]', 'Prec. Anomalia MAM - La Niña[2011]',\
           'Prec. Anomalia JJA - La Niña[2010]', 'Prec. Anomalia SON - La Niña[2010]']

# Criar a figura com 1 linha e 2 colunas
fig, axes = plt.subplots(2, 2, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

# Ajustar para que axes seja uma lista plana para iterar facilmente
axes = axes.flatten()

# Loop para configurar e plotar cada subplot
for i, ax in enumerate(axes):
    # Adicionar linhas de costa, limites de países e rios
    ax.coastlines(resolution='110m', linewidth=1)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    
    # Exemplo de adição de rios
    rivers_50m = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '50m')
    ax.add_feature(rivers_50m, facecolor='None', edgecolor='b', linewidth=1)
    
    # Contorno da Bacia Amazônica
    fname = 'amazonas/Limite_Cuenca_Amazonas_geogpsperu_juansuyo.shp'
    adm1_shapes = list(shpreader.Reader(fname).geometries())
    ax.add_geometries(adm1_shapes, ccrs.PlateCarree(), edgecolor='yellow', facecolor='none', 
                      alpha=1, linewidth=1)
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
    norm = mcolors.TwoSlopeNorm(vmin=-100, vcenter=0, vmax=100)
    # Exemplo de plotagem do campo de magnitude da anomalia de precipitação para cada mês (ajustar conforme necessário)
    anomaly_dataset_nina[i].plot(ax=ax, cmap=cmocean.cm.tarn, norm=norm, transform=ccrs.PlateCarree(),
                            levels=np.arange(-100,101,10), vmin=-100, vmax=100, 
                            cbar_kwargs={'label':'[$mm.mes^{-1}$]', 'ticks':np.arange(-100, 101, 20), 'shrink':0.75})

    # Plotar a topografia
    colorb = np.arange(1500, 1550, 50)
    ax.contour(ds_topo_ba.lon, ds_topo_ba.lat, ds_topo_ba['z'], 
               levels=colorb, colors=['brown'], linewidths=0.5)

    # Título para cada subplot (mês)
    ax.set_title(titles[i], fontsize=14)

# Ajustar o layout para evitar sobreposi��o
plt.subplots_adjust(hspace=0.45, wspace=0.15) #, left=0.125, bottom=0.1, right=0.9, top=0.9)
plt.tight_layout()
plt.show()
fig.savefig('imerg_anom_La-Nina_2010-2011.png', format='png', dpi = 300)

