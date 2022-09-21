# Algorithm Input and Output Data Definition (IODD)
[//]: # "This section should include the following:"
[//]: # "- Define and describe the content and format of input data required for the Level-2 ATBD."
[//]: # "- Define and describe the content and format of output data required for the Level2 ATBD (NetCDF CF-complaint containers are assumed)."
[//]: # "- Define and describe the content and format of auxiliary data required for the Level-2 ATBD."
[//]: # "- Define and describe the content and format of ancillary data required for the Level-2 ATBD."
[//]: # "- Define any other aspect considered relevant to the ATBD IODD."
The processing of L2 SST requires both TB observations and ancillary data. Additional information such as e.g. a land mask, a sea ice product and information about RFI is needed to remove contaminated PMW observations. The input for the SST retrieval algorithm consists of TB's at 1.4, 6.9, 10.6, 18.7 and 36.5 GHz (vertical and horizontal polarization) and additional satellite information as well as surface WS. The output is subskin SST and optionally surface WS, as described in this document.


## Input data
| Field | Description | Shape/Amount |
| --- | --- | --- |
| L1B TB L-band | L1B Brightness Temperatures at 1.4 GHz (V and H polarization)| Full swath or section of it (Nscans, Npos) |
| L1B TB C-band | L1B Brightness Temperatures at 6.9 GHz (V and H polarization)| Full swath or section of it (Nscans, Npos) |
| L1B TB X-band | L1B Brightness Temperatures at 10.6 GHz (V and H polarization)| Full swath or section of it (Nscans, Npos) |
| L1B TB K-band | L1B Brightness Temperatures at 18.7 GHz (V and H polarization)| Full swath or section of it (Nscans, Npos) |
| L1B TB Ka-band | L1B Brightness Temperatures at 36.5 GHz (V and H polarization)| Full swath or section of it (Nscans, Npos) |


## Output data
| Field | Description | Shape/Amount |
| --- | --- | --- |
| L2 SST | Sea Surface Temperature | Full swath or section of it (Nscans, Npos) |
| L2 WS | Wind Speed (optional) | Full swath or section of it (Nscans, Npos) |

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
