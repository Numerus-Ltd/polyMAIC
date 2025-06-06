# polyMAIC

## Overview
Matching-Adjusted Indirect Comparison (MAIC) is a technique used to account for the differences in baseline characteristics between treatment groups in indirect comparisons, especially when individual patient data are available for one treatment, but only aggregate-level data exist for the comparator. The PolyMAIC method, proposed by Alsop and Pont<sup>1</sup>, uses fourth order polynomial functions to estimate individual patient weights (IPW). 

## Features 
This package:
1. Contains the simulated data used in  Alsop and Pont<sup>1</sup>.
1. Estimates IPWs.
2. Provides matching algorithm diagnostics.
3. Produces a histogram of re-scaled IPWs.
   
## Usage
### Installation  
Install the development version using:
remotes::install_github("Numerus Ltd/Polymaic")

This package has not been uploaded to the CRAN repository. 

### Dependencies 
The polymaic() package is dependant upon the following packages: 
+ nloptr() 
+ Hmisc()
+ tictoc()
+ haven()
+ stringr()
+ ggplot2()

## Bibliography 
1. Alsop JC, Pont LO. Matching-adjusted indirect comparison via a polynomial-based non-linear optimization method. J Comp Eff Res. 2022 Jun;11(8):551-561. doi: 10.2217/cer-2021-0266. Epub 2022 May 4. PMID: 35506464.
2. Signorovitch J, Wu E, Yu A et al. Comparative effectiveness without head-to-head trials: a method for matching-adjusted indirect
comparisons applied to psoriasis treatment with adalimumab oretanercept. Pharmacoeconomics 28(10), 935–945 (2010).
3. Johnson S.G. SLSQP algorithm in NLopt. https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#slsqp.
4. Johnson S (2008). The NLopt nonlinear-optimization package. https://github.com/stevengj/nlopt.

## Citation 
Developers: J. Wilson, C. Reynish, L. Pont, J. Alsop.<br />
Initial Upload Date: 06-06-2025 <br />
For any enquires please contact polymaic@numerus.com 
