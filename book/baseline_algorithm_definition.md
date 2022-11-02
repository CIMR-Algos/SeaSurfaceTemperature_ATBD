# Baseline Algorithm Definition


## Forward Model
N/A.


## CIMR Level-1b re-sampling approach
The re-sampling approach for the {term}`CIMR` Level-1b data is to be investigated. First, simple methods, such as nearest neighbor and simple averaging are investigated, whereafter more advanced methods will be investigated.


## Algorithm Assumptions and Simplifications
In version 1 of the retrieval algorithm, {term}`AMSR2` data will be used as input instead of {term}`CIMR` data, as it is not yet available. Because of this, the channel combination used will not be the full {term}`CIMR` suite, as {term}`AMSR2` does not include the 1.4 GHz channels, but a so-called {term}`CIMR`-like combination, as investigated in {cite:t}`Nielsen2021`.

## Level-2 end to end algorithm functional flow diagram
The {term}`CIMR` {term}`SST` retrieval algorithm consists of an optional {term}`WS` retrieval algorithm followed by an {term}`SST` retrieval algorithm. A flow diagram of the algorithm is shown in {numref}`fig-algo-flow-diag`.
```{figure} figures/algo_flow_diagram/algo_flow_diagram.png
----
align: center
name: fig-algo-flow-diag
----
Set-up of the {term}`PMW` {term}`SST` retrieval algorithm using {term}`CIMR` orbital data and ancillary data as input. $\boldsymbol{\mathrm{A}}$, $\boldsymbol{\mathrm{B}}$ and $\boldsymbol{\mathrm{C}}$ refer to the regression coefficients in Equations {eq}`eq_ws_global`, {eq}`eq_ws_local` and {eq}`eq_sst_global`, respectively.
```

## Retrieval algorithm
The retrieval algorithm is a statistically-based algorithm for retrieving {term}`SST` given satellite {term}`TB`s and ancillary data, such as e.g. {term}`WS`. The main option is to use retrieved {term}`WS` from the integrated {term}`OEM` retrieval. Two additional options exist; the use of ancillary {term}`WS`, from e.g. reanalysis data, or retrieval of {term}`WS` using a 2-stage {term}`WS` retrieval algorithm.

### Mathematical description

#### WS retrieval algorithm
Optionally, a statistical retrieval algorithm can be used to retrieve {term}`WS` given satellite {term}`TB`s and ancillary data. The {term}`WS` retrieval algorithm is a two-stage algorithm based on the multiple linear regression model from {cite:t}`Alerskans2020`. In the first stage, an initial estimate of {term}`WS` is retrieved using a so called global regression-based retrieval algorithm, i.e. the model uses one set of regression coefficient for all retrievals. In the second stage, a final estimate of {term}`WS` is obtained through the use of localized algorithms, such that different sets of regression coefficients are used for the retrievals. Here, the retrieved {term}`WS` from the 1st-stage retrieval is used to bin the data and regression coefficients are obtained for a set of pre-defined {term}`WS` intervals. Hence, the 2nd-stage retrieval algorithm is trained to perform well over restricted {term}`WS` domains.

##### 1st-stage: Global retrieval algorithm
The first stage of the {term}`WS` retrieval algorithm uses a global regression model to obtain regression coefficients based on all training examples in the training dataset. The {term}`WS` retrieval algorithm is based on the {term}`NOAA` {term}`AMSR2` {term}`WS` retrieval algorithm {cite:p}`Chang2015` and uses {term}`TB`s and Earth incidence angle ($\theta_{EIA}$) to obtain initial retrieved {term}`WS`, $WS_a$,
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
Localized {term}`WS` retrieval algorithms are used in the second stage of the {term}`WS` retrieval algorithm. Here, the algorithms are defined for local {term}`WS` bins, such that separate regression coefficients are derived for each {term}`WS` interval, where $WS_{a}$ is used to select the correct {term}`WS` bin. In order to obtain robust regression coefficients, the minimum number of matchups required were set to 100. The localized algorithms are defined for {term}`WS` in the interval 0 to 20 ms$^{-1}$, with a bin size of 1 ms$^{-1}$, which gives a total of 20 localized algorithms. Like the 1st-stage retrieval, {term}`TB` and incidence angle are used in the localized {term}`WS` algorithm to obtain a final estimate of {term}`WS`, $WS_r$,
 ```{math}
  :label: eq_ws_local
  WS_{rk} = b_{0k} + \sum_{i=1}^{N_{ch}} (b_{1ik} t_{i} + b_{2ik} t^{2}_{i} ) + b_{3k} \theta,
 ```
where $t_{i}$ is given by Equation {eq}`eq_t11` and $\theta$ is defined by Equation {eq}`eq_theta`. The index $i$ refers to the summation over all $N_{ch}$ channels included in the retrieval algorithm, $k$ refers to the reference {term}`WS` bin and the coefficients $b_{0}$ to $b_{3}$ are regression coefficients, here together referred to as $B$, determined using the least-squares method. The final retrieved {term}`WS` is found by performing a linear interpolation between $WS_{rk}$ and the {term}`WS` retrieved using the closest neighboring {term}`WS` algorithm in order to avoid discontinuities

```{math}
 :label: eq_ws_interpol
 WS_{r} = \sum_{k=k_{0}}^{k_{0}+1} w_{k-k_{0}} WS_{rk}
```
where the interpolation weights $w$ are given by $w_{0}=1-\alpha$ and $w_{1}=\alpha$, with $\alpha=\frac{WS_{\alpha}}{\Delta k}-k_{0}$ and $k_{0}=floor(\frac{WS_{r}}{\Delta k})$, where $\Delta k=1~ms^{-1}$ is the {term}`WS` bin size.


#### SST retrieval algorithm
A statistical retrieval algorithm is used to retrieve {term}`SST` given satellite {term}`TB`s and ancillary data. The {term}`SST` retrieval algorithm is based on multiple linear regression and uses a global regression model to obtain regression coefficients based on all training examples in the training dataset. The algorithm is inspired by the {term}`RSS` {term}`AMSR-E` {term}`SST` retrieval algorithm {cite}`Wentz2000b` and the retrieval algorithm from {cite:t}`Alerskans2020`, in which the {term}`SST` is retrieved using {term}`TB`s, Earth incidence angle, retrieved {term}`WS` and the relative angle between wind direction and satellite azimuth angle ($\phi_{rel}$)
 ```{math}
  :label: eq_sst_global
  SST_{r} = c_{0} + \sum_{i=1}^{N_{ch}} (c_{1i} t_{i} + c_{2i} t^{2}_{i} ) + c_{3} \theta + c_{4} WS_{r} + \sum_{j=1}^{2} [c_{5j} \cos (j \phi_{rel}) + c_{6j} \sin (j \phi_{rel}) ],
 ```
 where again $t_{i}$ is given by Equation {eq}`eq_t11` and $\theta$ is defined by Equation {eq}`eq_theta`. The index $i$ refers to the summation over all $N_{ch}$ channels included in the retrieval algorithm; 6.9, 10.6, 18.7 and 36.5 GHz (dual polarization), and the coefficients $c_{0}$ to $c_{6}$ are regression coefficients, here together referred to as $C$, determined using the least-squares method.


### Input data
In the initial phase, the {term}`ESA` {term}`CCI` {term}`MMD` will be used for algorithm development and tuning. The {term}`ESA` {term}`CCI` {term}`MMD`, previously described in {cite:t}`Nielsen2018` and {cite:t}`Alerskans2020`, contains {term}`TB`s from the {term}`AMSR-E` level 2A and {term}`AMSR2` level 1R swath data products {cite:p}`Ashcroft2013,Maeda2016`. It also includes quality controlled in situ {term}`SST` observations from the  International Comprehensive Ocean-Atmosphere DataSet version 2.5.1 {cite:p}`Woodruff2011` and the Met Office Hadley Centre Ensembles dataset version 4.2.0 {cite:p}`Good2013`. Additional data include reanalysis data from the ERA-Interim {cite:p}`Dee2011` and ERA5 reanalyses {cite:p}`Hersbach2020`. The {term}`MMD` includes temporally matched and collocated matchups from the period June 2002 - October 2011 and July 2012 - December 2016.

{term}`AMSR2` {term}`TB`s will be used for the initial algorithm development and validation. A {term}`CIMR`-like channel combination will be used, based on {cite:t}`Nielsen2021`. For this set-up, necessary {term}`PMW` observations are:
- {term}`AMSR2` {term}`TB` channels: 6.9, 10.6, 18.7 and 36.5 GHz (dual polarization)

Furthermore, the use of simulated {term}`CIMR` {term}`TB`s in the development of the algorithm will be investigated. This has the advantage of including the L-band {term}`TB`s, which are not available from {term}`AMSR2`.

The next phase will see the use of {term}`CIMR` orbital data and here the following input data are needed:
- {term}`CIMR` {term}`TB` channels: 1.4, 6.9, 10.6, 18.7 and 36.5 GHz (dual polarization)


### Output data
The outputs from the regression retrieval algorithm are:
- Sea surface temperature ($SST_{r}$), in Kelvin.
- Wind speed ($WS_{r}$), in $~ms^{1}$ (optional).

### Auxiliary data
Data that is used as a complement to the retrieval, such as for flagging:
- Sea ice product
- Distance to coast
- Sun glint information
- Information for {term}`RFI` flagging

### Ancillary data
Data necessary for the retrieval:
- Earth incidence angle (satellite zenith angle)
- Satellite azimuth angle
- Reanalysis surface winds

### Validation process
Validation of the Level-2 {term}`SST` product is based on comparison of retrieved {term}`SST` with collocated and temporally matched in situ observations. In the initial phase, the {term}`ESA` {term}`CCI` {term}`MMD` is used for evaluation of the {term}`CIMR` {term}`SST` algorithm performance. The metrics used for the validation are standard verification metrics such as bias, standard deviation, root mean square error (RMSE) and coefficient of correlation (r).

Furthermore, the validation process will also include validation on Picasso scenes.
