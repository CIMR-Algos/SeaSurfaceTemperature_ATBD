# Baseline Algorithm Definition
[//]: # "This section should include the following:"
[//]: # "- Describe the retrieval method to be used for the Level-2 product. The ATBD shall provide an explanation of the scientific background and of the basic logical processing model."
[//]: # "- Define and describe the Forward Models to be used including models that are able to represent the target surface and its dynamics in a realistic fashion."
[//]: # "- Define and justify the re-gridding/re-sampling approach that shall be used to pre-process CIMR L1b data prior to applying the specific Level-2 retrieval algorithm. The regridding/re-sampling approach will be different depending on the Level-2 product. (NB. Eventually, this information will be used in a generic CIMR re-gridding/re-sampling tool as part of a separate IPF implementation activity that includes a Level-2 pre-processing re-gridding/re-sampling function)."
[//]: # "- Review and define rigorous error propagation methods allowing full characterisation and traceability of the products' uncertainty within the Level-2 products. The Contractor shall collaborate with the MACRAD study team."
[//]: # "- Review and define algorithm assumptions and simplifications that are used in the ATBD highlighting the anticipated impact of each on the final performance assessment. In addition, potential mitigation activities that could address these aspects shall be provided."
[//]: # "- Identify and review requirements for auxiliary data required to implement the algorithm and discuss their availability"
[//]: # "- Develop an end-to-end algorithm functional flow diagram for the Level-2 ATBD."
[//]: # "- Include any other aspect considered relevant to this part of the ATBD."
[//]: # "- For the Functional description of each algorithm step:"
[//]: # "    - Provide a mathematical description of the algorithm step."
[//]: # "    - Define interfaces to the algorithm step."
[//]: # "    - Define the dependencies for the algorithm step."
[//]: # "    - Define the input data for the algorithm step."
[//]: # "    - Define the output data for the algorithm step."
[//]: # "    - Define the configurable parameters for the algorithm step."
[//]: # "    - Define the auxiliary data for the algorithm step."
[//]: # "    - Define the ancillary data for the algorithm step."
[//]: # "    - Define the calibration process for the algorithm step, if necessary."
[//]: # "    - Define the validation process for the algorithm step."
[//]: # "    - Define the Calibration Data Set (ACDAT) required to calibrate the Level-2 algorithm step (if necessary)."
[//]: # "    - Define the Validation Data Set (AVDAT) required to validate the Level-2 algorithm step."
[//]: # "    - Define a Test Data Set (TDS) to verify the basic functioning of the entire end-to-end algorithm chain."
[//]: # "    - Define any other aspect considered relevant to the algorithm step"

[//]: # "(sec_qual_flags)="
[//]: # "## L2P Flags and Quality Levels"
[//]: # ""
[//]: # "The PMW SST CDR retrievals follow the GHRSST GDS 2.0 data specification {cite}`GHRSST2010` or L2P and each retrieval was assigned a quality level to denote the quality of the retrieval. The definition of quality levels, together with corresponding checks and thresholds, are shown in {numref}`tab_qual_levs`."
[//]: # ""
[//]: # "The quality of the individual SST retrievals is represented by a quality level, ranging from 0 to 5. Quality level 0 denotes the lowest quality indicator level, which is assigned if no data is retrieved. Quality level 1 is the lowest quality level for retrievals, indicating retrievals of bad quality which should not be used, whereas quality level 5 is the highest quality level, only given to retrievals with the best quality. Retrievals are assigned quality level 1 if the input data is of bad quality or if the retrieval is compromised, e.g. due to atmospheric and surface effects. The following criteria decide if a retrieval is of quality level 1:"
[//]: # ""
[//]: # "- AMSR-E scan quality or channel quality indicates bad satellite data."
[//]: # "- Any brightness temperature is outside the normal range (0K $\lt T_{B}\lt$320 K)."
[//]: # "- Sea ice contamination."
[//]: # "- Coastal contamination."
[//]: # "- Contamination due to RFI (masked according to {ref}`sec_rfi_filt`)."
[//]: # "- Rain contamination ($T_{B18V}\ge$240 K)."
[//]: # "- Sun glitter contamination ($\phi_{SGA}\lt$25$^{o}$)."
[//]: # "- Cases where the atmospheric contribution exceeds the information from the surface, i.e. if the difference between the horizontal and vertical polarisation brightness temperatures for channel 18-36 GHz is negative."
[//]: # "- The retrieved WS is outside the accepted range (0 $~m~s^{-1}\le WS_{r} \le$20$~m~s^{-1}$)."
[//]: # "- The retrieved SST is outside the accepted range (-2$^{o}$C$\le SST_{r} \le$ 35$^{o}$C)."
[//]: # "- The retrieved SST deviates with more than 10 $^{o}$C from a background SST."
[//]: # ""
[//]: # "Quality level 2, which denotes the worst-quality yet usable retrievals, is assigned to retrievals with a total uncertainty greater than 1. In addition, the proximity to sea ice and land is also used to determine if the retrieval is of quality level 2. If the distance to sea ice is less than 200 km or if the distance to land is less than 40 km, the retrieval is classified as being of quality level 2. Quality level 3 to 5 are determined based solely on the retrieved total SST uncertainty. If the SST uncertainty is in the range 0.5-1 K, the retrieval is assigned quality level 3 (low quality), if it is in the range 0.35-0.5 K, the retrieval is assigned quality level 4 (acceptable quality) and if the uncertainty is 0.35 K or smaller, the retrieval is assigned quality level 5 (best quality)."
[//]: # ""
[//]: # " ```{table} Quality levels with corresponding checks and thresholds."
[//]: # " :name: tab_qual_levs"
[//]: # " | #      |  **Quality Description** | **Checks & Thresholds** |   "
[//]: # " | ------ |   ---------------------- |  ---------------------  |"
[//]: # " | 0      |           No data        |                         |"
[//]: # " | 1      |          Bad data        | Quality controls & various atmospheric & surface effects |"
[//]: # " | 2      | Worst quality usable data| $\epsilon_{SST_{r}}\ge$1, Proximity to sea ice, Proximity to land |"
[//]: # " | 3      |        Low quality       | 0.5$\lt\epsilon_{SST_{r}}\lt$1 |"
[//]: # " | 4      |   Acceptable quality     | 0.35$\lt\epsilon_{SST_{r}}\le$0.5 |"
[//]: # " | 5      |        Best quality      | $\epsilon_{SST_{r}}\le$0.35       |"
[//]: # "```"
[//]: # ""


## Forward Model
Not applicable now.


## CIMR Level-1b re-sampling approach
Add text here.


## Algorithm Assumptions and Simplifications
In version 0 of the retrieval algorithm, AMSR2 data will be used as input instead of CIMR data, as it is not yet available. Because of this, the channel combination used will not be the full CIMR suite, as AMSR2 does not include the 1.4 GHz channels, but a so-called CIMR-like combination.

## Level-2 end to end algorithm functional flow diagram
The CIMR SST retrieval algorithm consists of a WS retrieval algorithm followed by an SST retrieval algorithm. A flow diagram of the algorithm is shown in {numref}`fig-algo-flow-diag`.
```{figure} figures/algo_flow_diagram/algo_flow_diagram.png
----
align: center
name: fig-algo-flow-diag
----
Set-up of the PMW SST retrieval algorithm using CIMR orbital data and ancillary data as input. $\boldsymbol{\mathrm{A}}$, $\boldsymbol{\mathrm{B}}$ and $\boldsymbol{\mathrm{C}}$ refer to the regression coefficients in Equations {eq}`eq_ws_global`, {eq}`eq_ws_local` and {eq}`eq_sst_global`, respectively.
```

## Retrieval algorithm
The retrieval algorithm is a statistically-based algorithm for retrieving SST given satellite brightness temperatures and ancillary data, such as e.g. wind speed. The main option is to use retrieved wind speed from the integrated optimal estimation model (OEM) retrieval. Two additional options exist; the use of ancillary wind speed, from e.g. reanalysis data, or retrieval of wind speed using a 2-stage wind speed retrieval algorithm.

### Mathematical description

#### WS retrieval algorithm
Optionally, a statistical retrieval algorithm can be used to retrieve WS given satellite brightness temperatures and ancillary data. The WS retrieval algorithm is the two-stage algorithm based on multiple linear regression from {cite:t}`Alerskans2020`. In the first stage, an initial estimate of WS is retrieved using a so called global regression-based retrieval algorithm, i.e. the model uses one set of regression coefficient for all retrievals. In the second stage, a final estimate of WS is obtained through the use of localized algorithms, such that different sets of regression coefficients are used for the retrievals. Here, the retrieved wind speed from the 1st-stage retrieval is used to bin the data and regression coefficients are obtained for a set of pre-defined wind speed intervals. Hence, the 2nd-stage retrieval algorithm is trained to perform well over restricted WS domains.

##### 1st-stage: Global retrieval algorithm
The first stage of the WS retrieval algorithm uses a global regression model to obtain regression coefficients based on all training examples in the training dataset. The WS retrieval algorithm is based on the NOAA AMSR2 WS retrieval algorithm {cite:p}`Chang2015` and uses brightness temperatures ($T_{B}$) and Earth incidence angle ($\theta_{EIA}$) to obtain initial retrieved WS, $WS_a$,
 ```{math}
  :label: eq_ws_global
  WS_{a} = a_{0} + \sum_{i=1}^{N_{ch}} (a_{1i} t_{i} + a_{2i} t^{2}_{i} ) + a_{3} \theta,
 ```
where
 ```{math}
 :label: eq_t11
 t_{i} = T_{Bi}-150
 ```
and
 ```{math}
 :label: eq_theta
 \theta = \theta_{EIA} - 55.
 ```

The summation index $i$ represents the summation over the $N_{ch}=8$ channels included in the retrieval algorithm: 6.9, 10.6, 18.7, and 36.5 (dual polarizations), and the coefficients $a_0$ to $a_3$ are the regression coefficients, here together referred to as $A$, determined using the least-squares method.

##### 2nd-stage: WS loclised retrieval algorithm
Localized WS retrieval algorithms are used in the second stage of the WS retrieval algorithm. Here, the algorithms are defined for local WS bins, such that separate regression coefficients are derived for each WS interval, where $WS_{a}$ is used to select the correct WS bin. In order to obtain robust regression coefficients, the minimum number of matchups required were set to 100. The localized algorithms are defined for WS in the interval 0 to 20 ms$^{-1}$, with a bin size of 1 ms$^{-1}$, which gives a total of 20 localized algorithms. Like the 1st-stage retrieval, brightness temperature and incidence angle are used in the localized WS algorithm to obtain a final estimate of WS, $WS_r$,
 ```{math}
  :label: eq_ws_local
  WS_{rk} = b_{0k} + \sum_{i=1}^{N_{ch}} (b_{1ik} t_{i} + b_{2ik} t^{2}_{i} ) + b_{3k} \theta,
 ```
where $t_{i}$ is given by Equation {eq}`eq_t11` and $\theta$ is defined by Equation {eq}`eq_theta`. The index $i$ refers to the summation over all $N_{ch}$ channels included in the retrieval algorithm, $k$ refers to the reference WS bin and the coefficients $b_{0}$ to $b_{3}$ are regression coefficients, here together referred to as $B$, determined using the least-squares method. The final retrieved WS is found by performing a linear interpolation between $WS_{rk}$ and the WS retrieved using the closest neighboring WS algorithm in order to avoid discontinuities

```{math}
 :label: eq_ws_interpol
 WS_{r} = \sum_{k=k_{0}}^{k_{0}+1} w_{k-k_{0}} WS_{rk}
```
where the interpolation weights $w$ are given by $w_{0}=1-\alpha$ and $w_{1}=\alpha$, with $\alpha=\frac{WS_{\alpha}}{\Delta k}-k_{0}$ and $k_{0}=floor(\frac{WS_{r}}{\Delta k})$, where $\Delta k=1~ms^{-1}$ is the WS bin size.


#### SST retrieval algorithm
A statistical retrieval algorithm is used to retrieve SST given satellite brightness temperatures and ancillary data. The SST retrieval algorithm is based on multiple linear regression and uses a global regression model to obtain regression coefficients based on all training examples in the training dataset. The algorithm is inspired by the Remote Sensing Systems (RSS) AMSR-E SST retrieval algorithm {cite}`Wentz2000b` and the retrieval algorithm from {cite:t}`Alerskans2020`, in which the SST is retrieved using brightness temperatures, Earth incidence angle, retrieved wind speed and the relative angle between wind direction and satellite azimuth angle ($\phi_{rel}$)
 ```{math}
  :label: eq_sst_global
  SST_{r} = c_{0} + \sum_{i=1}^{N_{ch}} (c_{1i} t_{i} + c_{2i} t^{2}_{i} ) + c_{3} \theta + c_{4} WS_{r} + \sum_{j=1}^{2} [c_{5j} \cos (j \phi_{rel}) + c_{6j} \sin (j \phi_{rel}) ],
 ```
 where again $t_{i}$ is given by Equation {eq}`eq_t11` and $\theta$ is defined by Equation {eq}`eq_theta`. The index $i$ refers to the summation over all $N_{ch}$ channels included in the retrieval algorithm; 6.9, 10.6, 18.7 and 36.5 GHz (dual polarization), and the coefficients $c_{0}$ to $c_{6}$ are regression coefficients, here together referred to as $C$, determined using the least-squares method.


[//]: # "### Quality flags"
[//]: # "Add later on?"

### Input data
In the initial phase, the ESA Climate Change Initiative (CCI) Multisensor Matchup Dataset (MMD) will be used for algorithm development and tuning. The ESA CCI MMD, previously described in {cite:t}`Nielsen2018` and {cite:t}`Alerskans2020`, contains brightness temperatures from the AMSR-E level 2A (L2A) and AMSR2 level 1R (L1R) swath data products {cite:p}`Ashcroft2013,Maeda2016`. It also includes quality controlled in situ SST observations from the International Comprehensive Ocean-Atmosphere DataSet (ICOADS) version 2.5.1 {cite:p}`Woodruff2011` and the Met Office Hadley Centre (MOHC) Ensembles dataset version 4.2.0 (EN4) {cite:p}`Good2013`. Additional data include reanalysis data from the ERA-Interim {cite:p}`Dee2011` and ERA5 reanalyses {cite:p}`Hersbach2020`. The MMD includes temporally matched and collocated matchups from the period June 2002 - October 2011 and July 2012 - December 2016.

AMSR2 brightness temperatures will be used for the initial algorithm development and validation. A CIMR-like channel combination will be used, based on {cite:t}`Nielsen2021`. For this set-up, necessary PMW observations are:
- AMSR2 TB channels: 6.9, 10.6, 18.7 and 36.5 GHz (dual polarization)

The next phase will see the use of CIMR orbital data and here the following input data are needed:
- CIMR TB channels: 1.4, 6.9, 10.6, 18.7 and 36.5 GHz (dual polarization)


### Output data
The outputs from the regression retrieval algorithm are:
- Sea surface temperature ($SST_{r}$), in Kelvin.
- Wind speed ($WS_{r}$), in $~ms^{1}$ (optional).

### Auxiliary data
Data that is used as a complement to the retrieval, such as for flagging:
- Sea ice product
- Distance to coast
- Information for RFI flagging

### Ancillary data
Data necessary for the retrieval:
- Earth incidence angle (satellite zenith angle)
- Satellite azimuth angle
- Reanalysis surface winds

### Validation process
Validation of the Level-2 SST product is based on comparison of retrieved SST with collocated and temporally matched in situ observations from drifting buoys. In the initial phase, the ESA CCI MMD is used for evaluation of the CIMR SST algorithm performance. The metrics used for the validation are standard verification metrics such as bias, standard deviation, root mean square error (RMSE) and coefficient of correlation (r).

Furthermore, the validation process will also include validation on Picasso scenes.
