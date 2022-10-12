# Level-2 product definition
[//]: # "This section should include the following:"
[//]: # "- Review and define the gridding scheme and map projection for the Level-2 product"
[//]: # "- Review and define the contents of the Level-2 product in Common Data Language (CDL) notation assuming NetCDF objects and data."
[//]: # "- Include any other aspect considered relevant to the Level-2 product definition."

The Level-2b SST product is obtained through application of the SST retrieval algorithm, given the Level-1b brightness temperatures. The product is provided in the native grid of the CIMR instrument, where each swath contains 135 samples per scan, `n_scans`, multiplied by the number of scans for each orbit, `n_samples_earth`. The product is provided in netCDF format and contains the variables listed in {numref}`l2-product-variables` .

```{table} Level-2 SST product variables
:name: l2-product-variables

| Variable name | Unit | Dimensions |
| --- | --- | --- |
| sea_surface_temperature | Kelvin | n_scans * n_samples_earth |
| wind_speed | ms$^{-1}$ | n_scans * n_samples_earth |
```
