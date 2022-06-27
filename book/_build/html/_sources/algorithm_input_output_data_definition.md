# Algorithm Input and Output Data Definition (IODD)


### Input data

The retrieval algorithm is designed to be able to use two sources of input data; a Multi-sensor Matchup Dataset (MMD), which is used for tuning, development and validation of the retrieval algorithm, and orbital AMSR-E and AMSR2 data, which is used for producing the PMW SST climate data record.

### ESA-CCI Multi-sensor Matchup Dataset
For tuning and development of the retrieval algorithm, as well as for assessment and validation of the performance of the algorithm, Multi-sensor Matchup Datasets (MMDs), version MMD6c and MMD6b, was used as input. The MMDs were generated using the Multi-sensor Matchup System (MMS) software which was developed within the ESA-CCI SST project and the European Union’s Horizon 2020 research and innovation programme under grant agreement No 638822 (FIDUECO project) {cite}`Block2018`. The MMD6c consists of AMSR-E orbital data matched to in situ SST measurements and MMD6b is the corresponding matchup database for AMSR2. The in situ dataset contains quality controlled measurements from Global Tropical Moored Buoy Array (GTMBA) data from NOAA PMEL {cite}`McPhaden2010`, the International Comprehensive Ocean-Atmosphere Dataset (ICOADS) version 2.5.1 {cite}`Woodruff2011` and the Met Office Hadley Centre (MOHC) Ensembles dataset version 4.2.0 (EN4) {cite}`Good2013`. In addition, the MMD also includes NWP data from ERA-Interim {cite}`Dee2011`, which have been interpolated in both space and time to the matchup location. Furthermore, the Cross-Calibrated Multi-Platform (CCMP) surface vector winds {cite}`Atlas2011` were collocated with the MMDs and was used for tuning and development, as well as validation, of the WS retrieval algorithm.

To obtain independent results, the MMDs are divided into seven subsets each, to be used for either tuning and development or validation of the algorithm:
- WS1_TRAIN: training subset for the 1st-stage WS retrieval algorithm
- WS1_TEST: validation subset for the 1st-stage WS retrieval algorithm
- WS2_TRAIN: training subset for the 2nd-stage WS retrieval algorithm
- WS2_TEST/SST_TRAIN: validation subset for the 2nd-stage WS retrieval algorithm, also used as training subset for the SST retrieval algorithm
- SST_TEST: validation subset for the SST retrieval algorithm
- UNCERT_TRAIN: training subset for the SST uncertainty retrieval algorithm
- UNCERT_TEST: validation subset for the SST uncertainty retrieval algorithm.

### AMSR-E orbital data
For producing a climate data record of PMW SST, the spatially resampled L2A swath data product AMSR-E V12 (Ashcroft et al., 2013), produced by Remote Sensing Systems (RSS) and distributed by NASA’s National Snow and Ice Data Center (NSIDC, https://nsidc.org/data/ae_l2a), is used as input. Here we used the brightness temperatures re-sampled to the 6.9 GHz resolutions. Hence the observations have a resolution footprint of 75 x 43 km, however, the data is distributed as a dataset with a spatial grid resolution of 10 km. Auxiliary data include NWP data from ERA-Interim. The Generalized Bayesian Cloud Screening (GBCS) software package (Mackie et al., 2010) is used to interpolate the NWP data to the satellite raster.

### AMSR2 orbital data
The spatially resampled AMSR2 L1R version 2 swath data product: Dataset of Brightness Temperature Modified Using the Antenna Pattern Matching Technique (Maeda et al., 2016) is used for producing the PMW SST CDR. This product contains similar spatially resampled brightness temperatures to the AMSR-E dataset. NWP data from ERA-Interim are used as auxiliary data. As for the AMSR-E processing, the GBCS software package is used to interpolate the auxiliary data to the satellite raster.

### Preprocessing
Preprocessing of the input data is necessary before running the regression model. Different fields need to be calculated depending on the input data used. For AMSR-E and AMSR2, the following fields need to be calculated:
- Relative angle between satellite azimuth angle and wind direction ($\phi_{REL}$), which is calculated as $\phi_{REL} = \phi_{SAT} - \phi_{WD}$, where $\phi_{SAT}$ is the satellite azimuth angle, relative to north, and $\phi_{WD}$ is the wind direction, relative to north, that the wind is blowing toward.
- Sun glint angle ($\phi_{SGA}$), which is calculated as $\phi_{SGA} = arccos(sin(\theta_{SOL})\times sin(\theta_{SAT})\times cos(\phi_{REL}+180)+cos(\theta_{SOL})\times cos(\theta_{SAT}))$, where $\theta_{SOL}$ and $\theta_{SAT}$ are the solar zenith angle and satellite zenith angle, respectively, and $\phi_{REL}$ is the relative azimuth angle between solar azimuth angle ($\phi_{SOL}$) and satellite azimuth angle ($\phi_{SAT}$).


### Output data

The outputs from the regression retrieval algorithm are:
- Sea surface temperature (SST_{r,baseline}), in Kelvin.
- Local systematic uncertainty component (\epsilon_{local}), in Kelvin.
- Random uncertainty component(\epsilon{rnd}), in Kelvin.
- Global systematic uncertainty component (\epsilon{global}), in Kelvin.
- Total SST uncertainty (\epsilon{SSTr}), in Kelvin.
- Wind speed (WS_{r}), in ms-1.
- L2P flags.
- Quality levels.

### Auxiliary data

Subsection text

### Ancillary data

Subsection text

