import cdsapi

dataset = "reanalysis-era5-pressure-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "specific_humidity",
        "u_component_of_wind",
        "v_component_of_wind"
    ],
    "pressure_level": [
        "100", "200", "250",
        "300", "400", "500",
        "600", "700", "800",
        "850", "900", "925",
        "1000"
    ],
    "year": [
        "2021", "2022", "2023",
        "2024"
        ],
    "month": [
    	"01", "02","03","04",
    	"05", "06","07","08",
    	"09", "10","11","12"
    	],
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [15, -85, -60, -30]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()

