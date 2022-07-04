### Baseline Algorithm Definition

Definition


### Retrieval Method

(sec_ws_retal)=
### Wind Speed Retrieval Algorithm

A regression-based retrieval algorithm is used to retrieve WS given satellite brightness temperatures and NWP fields. The WS retrieval algorithm described here is a two-step multiple linear regression model. In the first step, an initial estimate of WS is retrieved using a global retrieval algorithm, i.e. one set of regression coefficients is used for all wind speeds. In the second stage, a final WS is retrieved using specialised WS algorithms, i.e. the algorithm is trained to perform well over restricted WS domains.

### Step 1. Global Retrieval Algorithm

An initial estimate of wind speed ($WS_{\alpha}$) is obtained through the use of a global regression model, where the regression coefficients are obtained through training on the WSS1_TRAIN subsets. The WS retrieval algorithm is inspired by the NOAA AMSR-2 WS retrieval algorithm {cite}`Chang2015` and expresses WS in terms of brightness temperature ($T_{B}$) and incidence angle ($\theta_{EIA}$), as in Equation {eq}`eq_ws_glob`

 ```{math}
  :label: eq_ws_glob
  WS_{\alpha} = a_{0} + \sum_{i=1}^{10} a_{1i}t_{i} + a_{2i}t_{i}^{2}) + a_{3}\theta
  ```
  with
  ```{math}
  :label: eq_t11
  t_{i} = T_{Bi}-150
  ```
  for all channels except the 23.6 GHz channels. While for the two 23.6 GHz channels

  ```{math}
  :label: eq_t12
  t_{i} = ln(290-T_{Bi})
  ```
  and finally,
  ```{math}
  :label: eq_theta
  \theta = \theta_{EIA} - 55
  ```
The coefficients $\alpha_{0}$ to $\alpha_{3}$ are regression coefficients, referred to as $B_{global}$ determined through use of a training dataset (see Appendices A1-2), the summation index $i$ represents the summation over 10 AMSR-E channels, i.e. 6.9, 10.7, 18.7, 23.6 and 36.5 GHz (dual polarisation), while $T_{B}$ denotes the brightness temperature for the $i$th channel.

### Step 2. Specialised WS retrieval algorithms

In the second step, a final wind speed ($WS_{r}$) is retrieved through the use of a specialised WS regression model, using the same retrieval algorithm as in the first step. Regression coefficients are derived through training on subsets of the WS2\_TRAIN subsets, defined for restricted reference WS intervals. The retrieved WS from the first stage ($WS_{\alpha}$) is used to determine the correct WS bin, from which regression coefficients are selected to perform the WS retrieval. The specialized algorithms are derived for reference wind speeds in the interval 0 to 20$~m~s^{-1}$, with a bin size of 1$~m~s^{-1}$, giving a total of 20 specialised WS algorithms of the form as in Equation {eq}`eq_ws_spec`

```{math}
 :label: eq_ws_spec
 WS_{rk} = b_{0k} + \sum_{i=1}^{10} b_{1ik}t_{i} + b_{2ik}t_{i}^{2}) + b_{3k}\theta
```
where $t_{i}$ as in Equation {eq}`eq_t11`, $t_{i}$ as in Equation {eq}`eq_t12` and $\theta$ as in Equation {eq}`eq_theta` and $k$ denotes the reference WS. The coefficients $b_{0}$ to $b_{3}$ are regression coefficients, referred to as $B_{WS}$ determined through use of a training dataset (see Appendices B1-2). The final retrieved WS is found by performing a linear interpolation between $WS_{rk}$ and the WS retrieved using the closest neighbouring WS algorithm

```{math}
 :label: eq_ws_close
 WS_{r} = \sum_{k=k_{0}}^{k_{0}+1} w_{k-k_{0}} WS_{rk}
```
where $w_{0}=1-\alpha$, $w_{1}=\alpha$, $\alpha=\frac{WS_{\alpha}}{\Delta k}-k_{0}$, $k_{0}=floor(\frac{WS_{r}}{\Delta k})$ and $\Delta k=1~m~s^{-1}$ is the WS bin size.

(sec_sst_retal)=
### SST Retrieval Algoritm

A regression model to retrieve SST is used, ingesting satellite brightness temperatures and NWP fields. The retrieval algorithm is inspired by the RSS AMSR-E SST retrieval algorithm {cite}`Wentz2000b` and expresses SST in terms of brightness temperature ($T_{B}$), incidence angle ($\theta_{EIA}$), retrieved wind speed ($WS_{r}$) and the relative angle between satellite azimuth angle and NWP wind direction ($\phi_{REL}$).

The SST retrieval algorithm uses 12 brightness temperature channels; 6.9, 10.7, 18.7, 23.8, 36.5 and 89.0 (vertical and horizontal polarisation). The SST retrieval algorithm described here is a two-step multiple linear regression model with specialised regression algorithms, i.e trained to perform well over specialised domains. In the first stage, the algorithm is trained to perform well over restricted latitude domains and for ascending and descending orbit, respectively, whereas in the second stage, the algorithm is trained to perform well over restricted SST and WS domains.

### Specialised latitude and ascending/descending retrieval algorithms
An initial estimate of SST ($SST_{\alpha}$) is obtained through the use of a specialised orbit and latitude regression retrieval algorithm. Regression coefficients are obtained through training on subsets of the SST\_TRAIN subsets, defined for restricted reference latitude intervals and for ascending and descending orbits. Latitude and ascending or descending orbit are used to determine the correct latitude and orbit bin, from which regression coefficients are selected to perform the SST retrieval. The specialised algorithms are derived for reference latitudes in the interval -90$^{o}$ to 90$^{o}$, with a bin size of 2$^{o}$, and ascending (0) or descending (1) orbit, giving a total of 182 specialised latitude and orbit algorithms, as in Equation {eq}`eq_sst_alm`

 ```{math}
  :label: eq_sst_alm
  SST_{alm} = c_{0lm} + \sum_{i=1}^{12}(c_{1ilm}t_{i} + c_{2ilm}t_{i}^{2}) + c_{3lm}\theta + c_{4lm}WS_{r} + \sum_{j=1}^{2}(c_{5jlm} cos j \phi_{REL} + c_{6jlm} sin j \phi_{REL} )
  ```
, where $l$ denotes the reference latitude and $m$ denotes the reference orbit, which ranges from 0 (descending) to 1 (ascending). The coefficients $c_{0}$ through $c_{6}$ are regression coefficients, referred to as **$B_{LAT,ORB}$** (see Appendices C1-2) and are determined through the use of the SST\_TRAIN subset. The initial estimate of SST is found through linear interpolation between $SST_{alm}$ and the SST retrieved using the closest neighbouring latitude and orbit algorithm

```{math}
 :label: eq_sst_a1
 SST_{\alpha} = \sum_{l=l_{0}}^{l_{0}+1}w_{l-l_{0}} SST_{alm}
 ```
, where
```{math}
 :label: eq_sst_a2
 w_{0} = 1-\alpha
 w_{1} = \alpha
 \alpha = \frac{\phi_{LAT}}{\Delta l} - l_{0}
 l_{0} = floor(\frac{\phi_{LAT}}{\Delta l})
 ```
and $\phi_{LAT}$ denotes latitude, $\Delta l=2^{o}$.

### Specialised SST and WS retrieval algorithms   

In the second stage, final SST ($SST_{r}$) is retrieved through the use of a specialised SST and WS regression model. Regression coefficients are obtained through training on subsets of the SST\_TRAIN subset, defined for restricted reference SST and WS intervals. Retrieved wind speed ($WS_{r}$) and the retrieved SST from the first stage ($SST_{\alpha}$) are used to determine the correct SST and WS bin, from which regression coefficients are selected to perform the SST retrieval. The specialised algorithms are derived for reference SSTs in the interval -2$^{o}$C to 34$^{o}$C, with a bin size of 2$^{o}$C, and reference WS in the interval 0 to 20 $~m~s^{-1}$, with a bin size of 2 $~m~s^{-1}$, giving a total of 209 specialised SST and WS algorithms, as in Equation {eq}`eq_sst_spec`

```{math}
 :label: eq_sst_spec
 SST_{rnp} = d_{0np} + \sum_{i=1}^{12}(d_{1inp}t_{i} + d_{2inp}t_{i}^{2}) + d_{3np}\theta + d_{4np}WS_{r} + \sum_{j=1}^{2}(d_{5jnp} cos j \phi_{REL} + d_{6jnp} sin j \phi_{REL} )
 ```
where $n$ denotes the reference SST and $p$ denotes the reference WS. The coefficients $d_{0}$ through $d_{6}$ are regression coefficients, referred to as **$B_{SST,WS}$**, and are determined through use of subsets of the SST\_TRAIN subsets (see Appendices D1-2). The final retrieved SST is found by performing a bi-linear interpolation between $SST_{rnp}$ and the SSTs retrieved using the three closest neighbouring SST and WS algorithms as in Equation {eq}`eq_sst_r1`

```{math}
 :label: eq_sst_r1
 SST_{r} = \sum_{n=n_{0}}^{n_{0}+1}\sum_{p=p_{0}}^{p_{0}+1} \omega_{n-n_{0},p-p_{0}} SST_{rnp}
 ```
, where
```{math}
 :label: eq_sst_r2
 \omega_{0,0} = (1-\beta)\times(1-\gamma)
 \omega_{1,0} = \beta \times(1-\gamma)
 \omega_{0,1} = (1-\beta)\times \gamma
 \omega_{1,1} = \beta \times \gamma
 \beta = \frac{SST_{\alpha}}{\Delta n})-n_{0}, \gamma = \frac{WS_{r}}{\Delta p})-p_{0}
 n_{0} = floor(\frac{SST_{\alpha}}{\Delta n}), p_{0} = floor(\frac{WS_{r}}{\Delta p})
 ```
and $\Delta n = 2^{o}C$, $\Delta p = 2 $~m~s^{-1}$ is the SST and WS bin size, respectively.

### RFI Filter

The baseline SST retrieval algorithm, described in {numref}`sec_sst_retal` uses 12 brightness temperature channels; 6.9, 10.7, 18.7, 23.6, 36.5 and 89.0 GHz (dual polarisation). For detection of RFI, two additional SST retrieval algorithms were defined; the -10GHz and -18GHz algorithms. They are formulated exactly as the baseline algorithm, with the exception that they use only 10 brightness temperature channels; the same as the baseline algorithm minus the 10 GHz channel (-10GHz algorithm) and minus the 18 GHz channels (-18GHz algorithm). As for the baseline retrieval algorithm, WS is first retrieved using the two-step regression model with the baseline WS algorithm and then the two-step regression model is used to retrieve SST. See Appendices E-H for regression coefficients **$B_{LAT,ORB}$** and **$B_{SST,WS}$** for the -10GHz and -18GHz algorithms.

A new RFI mask, based on the two additional retrieval algorithms, has been developed. A 3$\sigma$-filter on the difference between retrieved SST for the two additional algorithms, -10GHz and -18GHz, and the baseline algorithm is used to detect RFI. Data are flagged if any of the following expressions {eq}`eq_sst_bsl1` and {eq}`eq_sst_bsl2` is true.

```{math}
 :label: eq_sst_bsl1
 SST_{r,baseline} - SST_{r,-10GHz} - \mu_{-10GHz} > 3\sigma_{-10GHz}
```

 ```{math}
  :label: eq_sst_bsl2
  SST_{r,baseline} - SST_{r,-18GHz} - \mu_{-18GHz} > 3\sigma_{-18GHz}
 ```

Here $SST_{r,-10GHz}$, $SST_{r,-18GHz}$ and $SST_{r,baseline}$ are the final retrieved SST using the -10GHz, - 18GHz and baseline algorithms, respectively, while $\mu_{-10GHz}$ and $\mu_{-18GHz}$ denote the mean difference $SST_{r,baseline} - SST_{r,-10GHz}$ and $SST_{r,baseline} - SST_{r,-18GHz}$ correspondingly, whereas $\sigma_{-10GHz}$ and $\sigma_{-18GHz}$ are the standard deviations of the differences. The $\mu$ and $\sigma$ used are shown in {numref}`tab_mean_sigma`  

 ```{table}
 :name: tab_mean_sigma
 | Sensor | $\mu_{-10GHz}$ (K) | $\mu_{-18GHz}$ (K) | $\sigma_{-10GHz}$ (K) | $\sigma_{-18GHz}$ (K) |
 | ---------- | ------- | ------ | ----- | ----- |
 | **AMSR-E** |  0.0024 | 0.0071 | 0.192 | 0.138 |
 | **AMSR2**  | -0.0087 | 0.0043 | 0.170 | 0.130 |
 ```

### Uncertainty Model



### Forward Model

Subsection Text


### CIMR Level-1b re-sampling approach

Subsection Text


### Algorithm Assumptions and Simplifications

Subsection Text

### Level-2 end to end algorithm functional flow diagram

The flow diagram of the algoritm is shown in {numref}`fig-algo-flow-diag`

```{figure} figures/cimr_algo_flow_diag.png
----
name: fig-algo-flow-diag
----
Setup of the DMI regression model for PMW SST retrievals using AMSR-E/2 orbital data as input. i denotes the algorithm used to retrieve SST; baseline (i=0), -10GHz (i=1) and -18GHz (i=2).
```

### Functional description of each Algorithm step

Subsection Text

### Mathematical description

SubSubsection Text

### Input data

The retrieval algorithm is designed to be able to use two sources of input data; a Multi-sensor Matchup Dataset (MMD), which is used for tuning, development and validation of the retrieval algorithm, and orbital AMSR-E and AMSR2 data, which is used for producing the PMW SST climate data record.

### Output data

The outputs from the regression retrieval algorithm are:
- Sea surface temperature ($SST_{r,baseline}$), in Kelvin.
- Local systematic uncertainty component ($\epsilon_{local}$), in Kelvin.
- Random uncertainty component($\epsilon{rnd}$), in Kelvin.
- Global systematic uncertainty component ($\epsilon{global}$), in Kelvin.
- Total SST uncertainty ($\epsilon{SSTr}$), in Kelvin.
- Wind speed ($WS_{r}$), in $~m~s^{1}$.
- L2P flags.
- Quality levels.

### Auxiliary data

SubSubsection Text

### Ancillary data

SubSubsection Text

##### Validation process

SubSubsection Text
