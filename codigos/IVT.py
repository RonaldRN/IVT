'''
Código que faz o plot do Transporte de umidade verticalmente integrado (IVT),
Aplica as funções do IVT do file IVT_Functions.py
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
from IVT_Functions import IVT, IVT_clim_season, IVT_clim_monthly, IVT_monthly

# Carregar o arquivo .nc"
ds = xr.open_dataset('nino-nina_2021-2024.nc')
ds = ds.sel(latitude=slice(15,-40), longitude=slice(-85,-30))

#print('-----------------------------------------------------------------')
# Ler arquivos mensais .nc 2001-2020
ds_monthly_2001_2020 = xr.open_dataset('AS_era5-montly_2001-2020.nc')
ds_monthly_2001_2020 = ds_monthly_2001_2020.sel(latitude=slice(15,-40), longitude=slice(-85,-30))

# Climatologia
#ivt_dataset_monthly = IVT_clim_monthly(ds_monthly_2001_2020)
# Salvar arquivo .nc
#ivt_dataset_monthly.to_netcdf("ivt_climatology_monthly.nc") 

# Sazonal
#ivt_dataset_season = IVT_clim_season(ds_monthly_2001_2020)
# Salvar arquivo .nc
#ivt_dataset_season.to_netcdf("ivt_climatology_season.nc") 

# IVT mensal desde Janeiro 2001 até Dezembro 2020
#ivt_dataset_monthly = IVT_monthly(ds_monthly_2001_2020)
# Salvar arquivo .nc
#ivt_dataset_monthly.to_netcdf("ivt_monthly_2001-2020.nc") 

# Ler arquivo da climatologia sazonal
ds_clim_season = xr.open_dataset('ivt_climatology_season.nc')
ds_clim_season = ds_clim_season.sel(latitude=slice(15,-40), longitude=slice(-85,-30))

# Ler arquivo da climatologia mensal
ds_clim_monthly = xr.open_dataset('ivt_climatology_monthly.nc')
ds_clim_monthly = ds_clim_monthly.sel(latitude=slice(15,-40), longitude=slice(-85,-30))

# Dataset da topografia
ds_topo = xr.open_dataset('topo_25.1.nc')
# Fazendo um slice para a �rea da Bacia Amazônica
ds_topo_ba = ds_topo.sel(lat=slice(-40, 10), lon=slice(-85, -40))

# ------------------------------------------------------------------------------------------------
# Figura da Climatologia [2001-2020] Mensal
# ------------------------------------------------------------------------------------------------
titles = ['IVT Janeiro [2001-2020]', 'IVT Fevereiro [2001-2020]', 'IVT Março [2001-2020]',
          'IVT Abril [2001-2020]', 'IVT Maio [2001-2020]', 'IVT Junho [2001-2020]', 'IVT Julho [2001-2020]',
          'IVT Agosto [2001-2020]', 'IVT Setembro [2001-2020]', 'IVT Outubro [2001-2020]',
          'IVT Novembro [2001-2020]', 'IVT Dezembro [2001-2020]']

# Criar a figura com 1 linha e 2 colunas
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
    ds_clim_monthly['IVT'][i].plot(ax=ax, cmap=cmocean.cm.rain, transform=ccrs.PlateCarree(),
                                   levels=np.arange(0,451,50), vmax=450, cbar_kwargs={'label':'[$kg.m^{-1}s^{-1}$]', 'shrink':0.5})

    # Agregar os dados de vento horizontal no formato de Vetor
    vetor = ax.quiver(ds_clim_monthly.longitude, ds_clim_monthly.latitude, ds_clim_monthly['IVT_u'][i],\
        ds_clim_monthly['IVT_v'][i], width=0.004, headlength=6, transform= ccrs.PlateCarree(), 
        regrid_shape=30, color='black')

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
fig.savefig('ivt_clim.png', format='png', dpi = 300)

# -------------------------------------------------------------------------------------------------
# Figura da Climatologia [2001-2020] Sazonal 
# ------------------------------------------------------------------------------------------------

titles = ['IVT DJF [2001-2020]', 'IVT MAM [2001-2020]', 'IVT JJA [2001-2020]', 'IVT SON [2001-2020]' ]

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

    # Exemplo de plotagem do campo de magnitude do IVT para cada mês (ajustar conforme necessário)
    ds_clim_season['IVT'][i].plot(ax=ax, cmap=cmocean.cm.rain, transform=ccrs.PlateCarree(),
                                   levels=np.arange(0,451,50), vmax=450, cbar_kwargs={'label':'[$kg.m^{-1}s^{-1}$]', 'shrink':0.5})

    # Agregar os dados de vento horizontal no formato de Vetor
    vetor = ax.quiver(ds_clim_season.longitude, ds_clim_season.latitude, ds_clim_season['IVT_u'][i],\
        ds_clim_season['IVT_v'][i], width=0.004, headlength=6, transform= ccrs.PlateCarree(), 
        regrid_shape=30, color='black')

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
fig.savefig('ivt_season.png', format='png', dpi = 300)

# -------------------------------------------------------------------------------------------------
# Selecionar Anos El Niño e La Niña
# -------------------------------------------------------------------------------------------------

# El Niño 2015/2016
ivt_dataset_nino_DJF = IVT(ds_monthly_2001_2020, date1=20151201, date2=20160101, date3=20160201)
ivt_dataset_nino_MAM = IVT(ds_monthly_2001_2020, date1=20160301, date2=20160401, date3=20160501)
ivt_dataset_nino_JJA = IVT(ds_monthly_2001_2020, date1=20150601, date2=20150701, date3=20150801)
ivt_dataset_nino_SON = IVT(ds_monthly_2001_2020, date1=20150901, date2=20151001, date3=20151101)

# La Niña 2010/2011
ivt_dataset_nina_DJF = IVT(ds_monthly_2001_2020, date1=20101201, date2=20110101, date3=20110201)
ivt_dataset_nina_MAM = IVT(ds_monthly_2001_2020, date1=20110301, date2=20110401, date3=20110501)
ivt_dataset_nina_JJA = IVT(ds_monthly_2001_2020, date1=20100601, date2=20100701, date3=20100801)
ivt_dataset_nina_SON = IVT(ds_monthly_2001_2020, date1=20100901, date2=20101001, date3=20101101)

# -------------------------------------------------------------------------------------------------
# Figura do IVT El Niño 2015/2016
# ------------------------------------------------------------------------------------------------

titles = ['IVT DJF - El Niño[2015-2016]', 'IVT MAM - El Niño[2016]', 'IVT JJA - El Niño[2015]', 'IVT SON - El Niño[2015]']
ivt_dataset_nino = [ivt_dataset_nino_DJF['IVT'], ivt_dataset_nino_MAM['IVT'], ivt_dataset_nino_JJA['IVT'], ivt_dataset_nino_SON['IVT']]
ivt_u_nino = [ivt_dataset_nino_DJF['IVT_u'], ivt_dataset_nino_MAM['IVT_u'], ivt_dataset_nino_JJA['IVT_u'], ivt_dataset_nino_SON['IVT_u']]
ivt_v_nino = [ivt_dataset_nino_DJF['IVT_v'], ivt_dataset_nino_MAM['IVT_v'], ivt_dataset_nino_JJA['IVT_v'], ivt_dataset_nino_SON['IVT_v']]

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

    # Exemplo de plotagem do campo de magnitude do IVT para cada mês (ajustar conforme necessário)
    # Plotar o campo de magnitude do IVT
    ivt_dataset_nino[i].plot(ax=ax, cmap=cmocean.cm.rain, transform=ccrs.PlateCarree(),
                        levels=np.arange(0,451,50), vmax=450, cbar_kwargs = {'label':'[$kg.m^{-1}s^{-1}$]','shrink':0.50})

    # Agregar os dados de vento horizontal no formato de Vetor
    vetor = ax.quiver(ivt_dataset_nino_DJF.longitude, ivt_dataset_nino_DJF.latitude, ivt_u_nino[i],\
        ivt_v_nino[i], width=0.004, headlength=6, transform= ccrs.PlateCarree(), 
        regrid_shape=30, color='black')

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
fig.savefig('ivt_nino.png', format='png', dpi = 300)

# -------------------------------------------------------------------------------------------------
# Figura do IVT La Niña 2010/2011
# ------------------------------------------------------------------------------------------------

titles = ['IVT DJF - La Niña[2010-2011]', 'IVT MAM - La Niña[2011]', 'IVT JJA - La Niña[2010]', 'IVT SON - La Niña[2010]']
ivt_dataset_nina = [ivt_dataset_nina_DJF['IVT'], ivt_dataset_nina_MAM['IVT'], ivt_dataset_nina_JJA['IVT'], ivt_dataset_nina_SON['IVT']]
ivt_u_nina = [ivt_dataset_nina_DJF['IVT_u'], ivt_dataset_nina_MAM['IVT_u'], ivt_dataset_nina_JJA['IVT_u'], ivt_dataset_nina_SON['IVT_u']]
ivt_v_nina = [ivt_dataset_nina_DJF['IVT_v'], ivt_dataset_nina_MAM['IVT_v'], ivt_dataset_nina_JJA['IVT_v'], ivt_dataset_nina_SON['IVT_v']]

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

    # Exemplo de plotagem do campo de magnitude do IVT para cada mês (ajustar conforme necessário)
    # Plotar o campo de magnitude do IVT
    ivt_dataset_nina[i].plot(ax=ax, cmap=cmocean.cm.rain, transform=ccrs.PlateCarree(),
                        levels=np.arange(0,451,50), vmax=450, cbar_kwargs = {'label':'[$kg.m^{-1}s^{-1}$]','shrink':0.50})

    # Agregar os dados de vento horizontal no formato de Vetor
    vetor = ax.quiver(ivt_dataset_nina_DJF.longitude, ivt_dataset_nina_DJF.latitude, ivt_u_nina[i],\
        ivt_v_nina[i], width=0.004, headlength=6, transform= ccrs.PlateCarree(), 
        regrid_shape=30, color='black')

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
fig.savefig('ivt_nina.png', format='png', dpi = 300)

# -------------------------------------------------------------------------
# Anomalias El Niño
# -------------------------------------------------------------------------
titles = ['IVT Anomalia DJF - El Niño[2015-2016]', 'IVT Anomalia MAM - El Niño[2016]',\
           'IVT Anomalia JJA - El Niño[2015]', 'IVT Anomalia SON - El Niño[2015]']

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

    # Exemplo de plotagem do campo de magnitude do IVT para cada mês (ajustar conforme necessário)
    (ivt_dataset_nino[i] - ds_clim_season['IVT'][i]).plot(ax=ax, cmap=cmocean.cm.diff_r, transform=ccrs.PlateCarree(),
                                levels=np.arange(-100,101,10), vmin=-100, vmax=100,
                                cbar_kwargs={'label':'[$kg.m^{-1}s^{-1}$]', 'shrink':0.5})

    # Agregar os dados de vento horizontal no formato de Vetor
    vetor = ax.quiver(ivt_dataset_nino_DJF.longitude, ivt_dataset_nino_DJF.latitude, (ivt_u_nino[i] - ds_clim_season['IVT_u'][i]),\
        (ivt_v_nino[i] - ds_clim_season['IVT_v'][i]), width=0.004, headlength=6, transform= ccrs.PlateCarree(), 
        regrid_shape=30, color='black')

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
fig.savefig('ivt_anom_El-Nino.png', format='png', dpi = 300)

# -------------------------------------------------------------------------
# Anomalias La Niña
# -------------------------------------------------------------------------
titles = ['IVT Anomalia DJF - La Niña[2010-2011]', 'IVT Anomalia MAM - La Niña[2011]',\
           'IVT Anomalia JJA - La Niña[2010]', 'IVT Anomalia SON - La Niña[2010]']

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

    # Exemplo de plotagem do campo de magnitude do IVT para cada mês (ajustar conforme necessário)
    (ivt_dataset_nina[i] - ds_clim_season['IVT'][i]).plot(ax=ax, cmap=cmocean.cm.diff_r, transform=ccrs.PlateCarree(),
                                levels=np.arange(-100,101,10), vmin=-100, vmax=100, 
                                cbar_kwargs={'label':'[$kg.m^{-1}s^{-1}$]', 'shrink':0.5})

    # Agregar os dados de vento horizontal no formato de Vetor
    vetor = ax.quiver(ivt_dataset_nina_DJF.longitude, ivt_dataset_nina_DJF.latitude, (ivt_u_nina[i] - ds_clim_season['IVT_u'][i]),\
        (ivt_v_nina[i] - ds_clim_season['IVT_v'][i]), width=0.004, headlength=6, transform= ccrs.PlateCarree(), 
        regrid_shape=30, color='black')

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
fig.savefig('ivt_anom_La-Nina.png', format='png', dpi = 300)
