# Introduction, purpose and scope

The sea surface temperature (SST) is an essential climate variable that is fundamental for climate monitoring and for the understanding of the air-sea interactions.
It has been observed from infrared satellites since the early 1980’s, but these observations are limited by clouds and severely impacted by aerosols {cite}`Reynolds1993`,{cite}`Reynolds2002`,{cite}`Vasquez2004`.
SST observations from passive microwave (PMW) observations are widely recognised as an important alternative to the infrared observations {cite}`Donlon2007`,{cite}`Donlon2009`.

They are not limited by clouds and the impact of aerosols is small {cite}`Chelton2005`,{cite}`Wentz2000a`.
This algorithm theoretical basis document (ATBD) describes in detail the DMI regression algorithm for the retrieval of SST from JAXA’s Advanced Microwave Scanning Radiometer - Earth Observing System (AMSR-E) and its follow-on instrument AMSR2 (Advanced Microwave Scanning Radiometer 2).

The algorithm has been used within the European Space Agency Climate Change Initiative Sea Surface Temperature (ESA-CCI SST) project to generate a global climate data record (CDR) of level 2 SSTs with associated uncertainties.
In addition, the retrieval also includes wind speeds (WSs). A consistent algorithm has been used for both the AMSR-E and AMSR2 observations.
