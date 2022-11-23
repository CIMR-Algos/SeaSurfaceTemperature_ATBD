# Level-2 product definition
The Level-2b {term}`SST` product is obtained through application of the {term}`SST` retrieval algorithm, given the Level-1b {term}`TB`s. The product is provided in the native grid of the {term}`CIMR` instrument, where each swath contains `n_samples` samples per scan, multiplied by the number of scans for each orbit, `n_samples_earth`. The product is provided in netCDF format and contains the variables listed in {numref}`l2-product-variables` .

```{table} Level-2 SST product variables
:name: l2-product-variables

| name | description | units | dimensions |
| --- | --- | --- | --- |
| sea_surface_temperature | sea surface temperature | Kelvin | n_scans * n_samples_earth |
| sea_surface_temperature uncertainty | sea surface temperature uncertainty | Kelvin | n_scans * n_samples_earth |
| status_flag | A flag indicating the quality of the SST retrieval | n/a | n_scans * n_samples_earth |
| wind_speed | surface wind speed | ms$^{-1}$ | n_scans * n_samples_earth |
```
