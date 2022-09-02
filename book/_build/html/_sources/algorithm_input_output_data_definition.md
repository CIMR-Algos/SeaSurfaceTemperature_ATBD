# Algorithm Input and Output Data Definition (IODD)
%This section should include the following:
%- Define and describe the content and format of input data required for the Level-2 ATBD.
%- Define and describe the content and format of output data required for the Level2 ATBD (NetCDF CF-complaint containers are assumed).
%- Define and describe the content and format of auxiliary data required for the Level-2 ATBD.
%- Define and describe the content and format of ancillary data required for the Level-2 ATBD.
%- Define any other aspect considered relevant to the ATBD IODD.


## Input data
- CIMR TB channels: 1.4, 6.9, 10.6, 18.7 and 36.5 GHz (dual polarization)
- Earth incidence angle (satellite zenith angle)
- Satellite azimuth angle
- NWP surface winds

### Preprocessing
The following pre-processing of input data is necessary:
- Calculation of NWP (meteorological) wind direction ($\phi_{WD}$) from the zonal ($u$) and meridional ($v$) surface wind components, following:
```{math}
:label: eq-preprocess-winddir
\phi_{WD} = 180 + \frac{180}{\pi} atan2(u,v),
```
- Calculation of relative angle ($\phi_{REL}$) between satellite azimuth angle ($\phi_{SAT}$) and NWP wind direction, following:
```{math}
:label: eq-preprocess-phirel
\phi_{REL} = \phi_{SAT} - \phi_{WD}
```

## Output data
- Retrieved SST
- Retrieved WS (optional)

## Auxiliary data
Data as a complement to the retrieval, e.g. for flagging:
- Sea ice product
- Distance to coast
- Information for RFI flagging

## Ancillary data
Data that is necessary for the retrieval:
- Earth incidence angle (satellite zenith angle)
- Satellite azimuth angle
- Reanalysis surface winds
