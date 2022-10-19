# Algorithm Input and Output Data Definition (IODD)
The processing of L2 {term}`SST` requires both {term}`TB` observations and ancillary data. Additional information such as e.g. a land mask, a sea ice product and information about {term}`RFI` is needed to remove contaminated {term}`PMW` observations. The input for the {term}`SST` retrieval algorithm consists of {term}`TB`s at 1.4, 6.9, 10.6, 18.7 and 36.5 GHz (vertical and horizontal polarization) and additional satellite information as well as surface {term}`WS`. The output is subskin {term}`SST` and optionally surface {term}`WS`, as described in this document.


## Input data
| Field | Description | Shape/Amount |
| --- | --- | --- |
| L1B TB L-band | L1B Brightness Temperatures at 1.4 GHz (V and H polarization)| Full swath or section of it (Nscans, Npos) |
| L1B TB C-band | L1B Brightness Temperatures at 6.9 GHz (V and H polarization)| Full swath or section of it (Nscans, Npos) |
| L1B TB X-band | L1B Brightness Temperatures at 10.6 GHz (V and H polarization)| Full swath or section of it (Nscans, Npos) |
| L1B TB Ku-band | L1B Brightness Temperatures at 18.7 GHz (V and H polarization)| Full swath or section of it (Nscans, Npos) |
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
- Sun glint information
- Information for {term}`RFI` flagging

## Ancillary data
Data that is necessary for the retrieval:
- Earth incidence angle (satellite zenith angle)
- Satellite azimuth angle
- Surface winds
