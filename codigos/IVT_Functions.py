# Funções para o cálculo do Transporte de umidade verticalmente integrado (IVT)

import numpy as np
import xarray as xr
import pandas as pd

# -------------------------------------------------------------------------------------------------
# Função que calcula a Climatologia sazonal [2001-2020] do 
# Integrated Verticallly Water Vapor Transport(IVT) 
# ------------------------------------------------------------------------------------------------

def IVT_clim_season(ds):
    '''
    Função IVT_clim_season cálcula a climatologia sazonal 
    do IVT (Integrated Vertically Vapor Water Transport)
    Inputs:
    - ds: dataset em formato NetCDF
    
    Outputs:
    - ivt_dataset: Dataset com o IVT, IVT_u, IVT_v
    '''

    ds = ds.sel(pressure_level=slice(1000,500))
    # Converter a coordenada 'date' para datetime64
    ds = ds.assign_coords(date=pd.to_datetime(ds['date'].astype(str), format='%Y%m%d'))
    # Calcular a climatologia de DJF [2001-2020]   
    ds_season = ds.groupby('date.season').mean(dim=['date'])

    # Obter os n�veis de press�o (normalmente em hPa)
    pressure = ds_season['pressure_level'].values * 100   # Converte hPa para Pa
    # Ajustando a dimens�o de pressure para (13, 1, 1)
    pressure = pressure[:, None, None]  # Adiciona duas dimens�es extras

    # Inicializar listas para armazenar resultados sazonais
    ivt_seasonal_list = []
    ivt_u_seasonal_list = []
    ivt_v_seasonal_list = []

    # Calcular IVT para cada estação
    for season in ['DJF', 'MAM', 'JJA', 'SON']:
        # Selecionar variáveis para a estação atual
        q_season = ds_season['q'].sel(season=season)
        u_season = ds_season['u'].sel(season=season)
        v_season = ds_season['v'].sel(season=season)
        
        # Constante de gravidade
        g = 9.81  # m/s²
        
        # Cálculo do IVT usando a regra dos trapézios
        IVT_u = -np.trapz(q_season * u_season, x=pressure, axis=0) / g
        IVT_v = -np.trapz(q_season * v_season, x=pressure, axis=0) / g
        IVT = np.sqrt(IVT_u**2 + IVT_v**2)
        
        # Adicionar resultados da estação à lista
        ivt_seasonal_list.append(IVT)
        ivt_u_seasonal_list.append(IVT_u)
        ivt_v_seasonal_list.append(IVT_v)

    # Cria um xarray.Dataset com os três arrays em um único dataset
    ivt_dataset = xr.Dataset({
        "IVT": (["season", "latitude", "longitude"], ivt_seasonal_list),
        "IVT_u": (["season", "latitude", "longitude"], ivt_u_seasonal_list),
        "IVT_v": (["season", "latitude", "longitude"], ivt_v_seasonal_list)
        }, 
        coords={
            "season": ['DJF', "MAM", "JJA", "SON"],
            "latitude": ds_season['latitude'],
            "longitude": ds_season['longitude']
        }
    )

    return ivt_dataset

# -------------------------------------------------------------------------------------------------
# Função que calculo a Climatologia mensal [2001-2020] do 
# Integrated Verticallly Water Vapor Transport(IVT) 
# ------------------------------------------------------------------------------------------------

def IVT_clim_monthly(ds):
    '''
    Função IVT_clim_monthly calcula a climatologia mensal 
    do IVT (Integrated Vertically Vapor Water Transport).
    Inputs:
    - ds: dataset em formato NetCDF com uma dimensão de tempo
    
    Outputs:
    - ivt_dataset: Dataset com o IVT, IVT_u, IVT_v
    '''
    ds = ds.sel(pressure_level=slice(1000,500))
    # Converter a coordenada 'date' para datetime64
    ds = ds.assign_coords(date=pd.to_datetime(ds['date'].astype(str), format='%Y%m%d'))
    # Calcular a climatologia de DJF [2001-2020]   
    ds_monthly = ds.groupby('date.month').mean(dim=['date'])
    # Obter os níveis de pressão (normalmente em hPa)
    pressure = ds_monthly['pressure_level'].values * 100   # Converte hPa para Pa
    pressure = pressure[:, None, None]  # Adiciona duas dimensões extras

    # Selecionar variáveis de umidade específica e componentes do vento
    q = ds_monthly['q']
    u = ds_monthly['u']
    v = ds_monthly['v']

    # Aceleração da gravidade (m/s²)
    g = 9.81

    # Calcular o IVT para cada mês
    ivt_monthly_list = []
    ivt_u_monthly_list = []
    ivt_v_monthly_list = []

    for month in range(1, 13):
        # Seleciona os dados do mês atual
        q_month = q.sel(month=int(month))#.mean(dim='time')
        u_month = u.sel(month=int(month))#.mean(dim='time')
        v_month = v.sel(month=int(month))#.mean(dim='time')

        # Calcula IVT para o mês usando a regra dos trapézios
        IVT_u = -np.trapz(q_month * u_month, x=pressure, axis=0) / g
        IVT_v = -np.trapz(q_month * v_month, x=pressure, axis=0) / g
        IVT = np.sqrt(IVT_u**2 + IVT_v**2)

        # Adiciona aos resultados mensais
        ivt_monthly_list.append(IVT)
        ivt_u_monthly_list.append(IVT_u)
        ivt_v_monthly_list.append(IVT_v)

    # Cria um xarray.Dataset com os três arrays em um único dataset
    ivt_dataset = xr.Dataset({
        "IVT": (["month", "latitude", "longitude"], ivt_monthly_list),
        "IVT_u": (["month", "latitude", "longitude"], ivt_u_monthly_list),
        "IVT_v": (["month", "latitude", "longitude"], ivt_v_monthly_list)
        }, 
        coords={
            "month": range(1, 13),
            "latitude": ds_monthly['latitude'],
            "longitude": ds_monthly['longitude']
        }
    )
    
    return ivt_dataset

# -------------------------------------------------------------------------------------------------
# Função que calculo o Integrated Verticallly Water Vapor Transport(IVT) 
# para períodos específicos
# ------------------------------------------------------------------------------------------------

# Funçãoo para períodos específicos
def IVT(ds, date1:int, date2:int, date3:int):

    '''
    Função IVT calcula o IVT (Integrated Vertically Vapor Water Transport).
    Inputs:
    - ds: dataset em formato NetCDF com uma dimensão de tempo
    - date1: Data 1
    - date2: Data 2
    - date3: Data 3
    
    Outputs:
    - ivt_dataset: Dataset com o IVT, IVT_u, IVT_v
    '''

    ds = ds.sel(pressure_level=slice(1000,500))
    # Obter os n�veis de press�o (normalmente em hPa)
    pressure = ds['pressure_level'].values * 100   # Converte hPa para Pa
    # Ajustando a dimens�o de pressure para (13, 1, 1)
    pressure = pressure[:, None, None]  # Adiciona duas dimens�es extras
    # Selecionar a vari�vel 'q' (q-Specific humidity [kg. kg-1])
    q_date1 = ds['q'].sel(date=date1)
    q_date2 = ds['q'].sel(date=date2)
    q_date3 = ds['q'].sel(date=date3)
    q_mean = (q_date1 + q_date2 + q_date3) / 3

    # Selecionar a vari�vel 'u' (U-component of wind [m.s-1])
    u_date1 = ds['u'].sel(date=date1)
    u_date2 = ds['u'].sel(date=date2)
    u_date3 = ds['u'].sel(date=date3)
    u_mean = (u_date1 + u_date2 + u_date3) / 3

    # Selecionar a vari�vel 'v' (V-component of wind [m.s-1])
    v_date1 = ds['v'].sel(date=date1)
    v_date2 = ds['v'].sel(date=date2)
    v_date3 = ds['v'].sel(date=date3)
    v_mean = (v_date1 + v_date2 + v_date3) / 3

    # Calcular o transporte integrado verticalmente (IVT)
    g = 9.81  # Acelera��o da gravidade (m/s�)

    # Integral aproximada usando a regra dos trap�zios
    IVT_u = -np.trapz(q_mean * u_mean, x=pressure, axis=0) / g
    IVT_v = -np.trapz(q_mean * v_mean, x=pressure, axis=0) / g

    # Vetor IVT
    IVT = np.sqrt(IVT_u**2 + IVT_v**2)

    # Converter para um xarray.DataArray com dimensões geogr�ficas (latitude, longitude)
    ivt_dataset = xr.Dataset(
        {
        "IVT": (["latitude", "longitude"], IVT),
        "IVT_u": (["latitude", "longitude"], IVT_u),
        "IVT_v": (["latitude", "longitude"], IVT_v),
        },
        coords={
            "latitude": ds['latitude'], 
            "longitude": ds['longitude']}
    )

    return ivt_dataset

# -------------------------------------------------------------------------------------------------
# Função que calculo o Integrated Verticallly Water Vapor Transport(IVT) 
# para cada tempo, longitude e latitude
# ------------------------------------------------------------------------------------------------

def IVT_monthly(ds):
    
    '''
    Função IVT_monthly calcula o IVT (Integrated Vertically Vapor Water Transport)
    para cada tempo do dataset (assumindo dados mensais).
    
    Inputs:
    - ds: dataset em formato NetCDF com dados médios mensais de 2001 até 2020
    
    Outputs:
    - ivt_dataset: Dataset com o IVT calculado para cada tempo, latitude e longitude
    '''

    ds = ds.sel(pressure_level=slice(1000,500))
    # Converter a coordenada 'date' para datetime64
    ds = ds.assign_coords(date=pd.to_datetime(ds['date'].astype(str), format='%Y%m%d'))

    # Obter os níveis de pressão (normalmente em hPa)
    pressure = ds['pressure_level'].values * 100  # Converte hPa para Pa
    # Ajustando a dimensão de pressure para (num_níveis, 1, 1) para facilitar broadcast
    pressure = pressure[:, None, None]  # Adiciona duas dimensões extras

    # Inicializar listas para armazenar resultados mensais
    ivt_list = []
    ivt_u_list = []
    ivt_v_list = []

    # Iterar sobre cada tempo no dataset
    for time_index in range(ds.dims['date']):
        # Selecionar variáveis para o tempo atual
        q_time = ds['q'].isel(date=time_index)
        u_time = ds['u'].isel(date=time_index)
        v_time = ds['v'].isel(date=time_index)
        
        # Constante de gravidade
        g = 9.81  # m/s²
        
        # Cálculo do IVT usando a regra dos trapézios ao longo da dimensão vertical (eixo 0)
        IVT_u = -np.trapz(q_time * u_time, x=pressure, axis=0) / g
        IVT_v = -np.trapz(q_time * v_time, x=pressure, axis=0) / g
        IVT = np.sqrt(IVT_u**2 + IVT_v**2)
        
        # Adicionar resultados do tempo atual às listas
        ivt_list.append(IVT)
        ivt_u_list.append(IVT_u)
        ivt_v_list.append(IVT_v)

    # Converter listas para arrays e empilhar ao longo da dimensão do tempo
    ivt_array = np.stack(ivt_list, axis=0)
    ivt_u_array = np.stack(ivt_u_list, axis=0)
    ivt_v_array = np.stack(ivt_v_list, axis=0)

    # Criar um xarray.Dataset com os três arrays em um único dataset
    ivt_dataset = xr.Dataset({
        "IVT": (["time", "latitude", "longitude"], ivt_array),
        "IVT_u": (["time", "latitude", "longitude"], ivt_u_array),
        "IVT_v": (["time", "latitude", "longitude"], ivt_v_array)
        }, 
        coords={
            "time": ds['date'].values,
            "latitude": ds['latitude'],
            "longitude": ds['longitude']
        }
    )

    return ivt_dataset