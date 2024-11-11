'''
Código que faz o plot das tendências do transporte de umidade verticalmente integrado (IVT)
calculado através do test estatístico de Mann-Kendall (previamente calculado)
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

# Ler arquivo das tendências mensais
ds_trend = xr.open_dataset('ivt_trend_2001-2020.nc')

# Dataset da topografia
ds_topo = xr.open_dataset('topo_25.1.nc')
# Fazendo um slice para a �rea da Bacia Amazônica
ds_topo_ba = ds_topo.sel(lat=slice(-40, 10), lon=slice(-85, -40))

# ------------------------------------------------------------------------------------------------
# Figura das tendências mensal Test Mann-Kendall [2001-2020] 
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
    #rivers_50m = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '50m')
    #ax.add_feature(rivers_50m, facecolor='None', edgecolor='b', linewidth=1)
    
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

    # Plot das tendências do IVT para cada mês (ajustar conforme necessário)
    ds_trend['slope'].isel(month=i).plot(ax=ax, cmap=cmocean.cm.balance_r, transform=ccrs.PlateCarree(),
                                   vmin=-6, vmax=6, cbar_kwargs={'label':'[$kg.m^{-1}s^{-1}$]', 'shrink':0.5, 'ticks': np.arange(-6, 7,2)})

    # Marcas as áreas onde o test é estatísticamente significativo ao nível de 95% com tendências positivas
    ax.contourf(ds_trend['longitude'], ds_trend['latitude'], ds_trend['numeric_trend'].isel(month=i) == 1,
                colors='none', hatches = ['XXXX'], levels = [0.5, 1.5])
    # Marcas as áreas onde o test é estatísticamente significativo ao nível de 95% com tendências negativas
    ax.contourf(ds_trend['longitude'], ds_trend['latitude'], ds_trend['numeric_trend'].isel(month=i) == -1, 
                colors='none', hatches = ['XXXX'], levels = [0.5, 1.5])
   
    # Plotar a topografia
    colorb = np.arange(1500, 1550, 50)
    ax.contour(ds_topo_ba.lon, ds_topo_ba.lat, ds_topo_ba['z'], 
               levels=colorb, colors=['brown'], linewidths=0.5)

    # Título para cada subplot (mês)
    ax.set_title(titles[i], fontsize=14)

# Ajustar o layout para evitar sobreposição
plt.subplots_adjust(hspace=0.5, wspace=0.15) #, left=0.125, bottom=0.1, right=0.9, top=0.9)
plt.tight_layout()
plt.show()
fig.savefig('ivt_trend.png', format='png', dpi = 300)