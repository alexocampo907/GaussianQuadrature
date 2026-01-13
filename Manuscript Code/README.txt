# README: Calculating True Values in Causal Inference Simulation Studies via Gaussian Quadrature

This repository contains the reproduction code for the manuscript: 
"Revealing the truth: Calculating True Values in Causal Inference Simulation Studies via Gaussian Quadrature"

---

## General Information

- **Software:** R version 4.4.3
- **Data Requirements:** None. All data are simulated within the scripts. Simulation results can alternatively be loaded under "Results/R Objects/".
- **Manual Edits Required:** None. The parameters load.sim=F and iter <- 10000 can be modified if desired.
- **Hardware Note:** Computation times estimated below are based on a standard workstation.

---

## File Structure & Contents

The following R scripts below reproduce the specific sections and figures cited in the paper. Each .R file can be run in isolation to generate both intermediate results and figures stored in the "Results" folder. 

No manual edits are required unless one would like toggle off the use the intermediate results to generate outputs and run the simulations. The load.sim=T loads the .rds files that contain the intermediated results in "Results/R Objects" can be swtiched to load.sim=F to rerun the original scripts. If one would like to rerun scripts with reduced runtime, the code line "iter <- 10000" in each of the scripts can be reduced. 

### 1. [3_Confounding.R]
- **Description:** Reproduces the results and figures presented in Section 3 regarding confounding bias and integration.
- **Output:** Reproduces Figures 4-8.
- **Estimated Runtime:** ~100 minutes. 
To reduce runtime, lower any text specifying the number of iterations (iter <- 10000) 


### 2. [4.1_CDE_continuous.R]
- **Description:** Performs Monte Carlo (MC) integration and Gauss-Hermite quadrature for the Controlled Direct Effect (CDE) with a continuous endpoint.
- **Output:** Reproduces the top panel of Figure 10.
- **Estimated Runtime:** ~11 minutes. 
To reduce runtime, lower any text specifying the number of iterations (iter <- 10000)

### 3. [4.1_CDE_binary.R]
- **Description:** Performs MC integration and Gauss-Hermite quadrature for the CDE with a binary endpoint.
- **Output:** Reproduces the bottom panel of Figure 10.
- **Estimated Runtime:** ~13 minutes. 
To reduce runtime, lower any text specifying the number of iterations (iter <- 10000)

### 4. [4.2_RMST.R]
- **Description:** Performs MC integration and Gauss-Hermite quadrature for mediation effects on the Restricted Mean Survival Time (RMST) scale.
- **Output:** Reproduces all panels of Figure 11.
- **Estimated Runtime:** ~8 minutes. 
To reduce runtime, lower any text specifying the number of iterations (iter <- 10000)

### 5. [Results/]
A local directory where all figures and numerical outputs are automatically saved upon execution of the scripts. Subdirectories are "Plots" for figures and "R Objects"" for storage of simulation results.

---

## Computational Environment & Dependencies

The following environment was used for the primary analysis. Key methodological packages include ` mvQuad` and `spatstat.random` for quadrature rules.

### Session Info:
Software details:

> sessionInfo()
R version 4.4.3 (2025-02-28)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gridExtra_2.3         spatstat.random_3.3-3 spatstat.geom_3.3-6   spatstat.univar_3.1-2
 [5] spatstat.data_3.1-6   survival_3.8-3        ggdist_3.3.2          gaussquad_1.0-3      
 [9] orthopolynom_1.0-6.1  mvQuad_1.0-8          pracma_2.4.4          lubridate_1.9.4      
[13] forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.4          
[17] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.2        
[21] tidyverse_2.0.0      

loaded via a namespace (and not attached):
 [1] generics_0.1.3       lattice_0.22-7       stringi_1.8.7        spatstat.utils_3.1-3
 [5] hms_1.1.3            magrittr_2.0.3       grid_4.4.3           timechange_0.3.0    
 [9] Matrix_1.7-3         scales_1.3.0         textshaping_1.0.0    cli_3.6.4           
[13] rlang_1.1.6          crayon_1.5.3         polyclip_1.10-7      munsell_0.5.1       
[17] splines_4.4.3        withr_3.0.2          tools_4.4.3          deldir_2.0-4        
[21] polynom_1.4-1        tzdb_0.5.0           colorspace_2.1-1     vctrs_0.6.5         
[25] R6_2.6.1             lifecycle_1.0.4      ragg_1.4.0           pkgconfig_2.0.3     
[29] pillar_1.10.2        gtable_0.3.6         data.table_1.17.0    glue_1.8.0          
[33] Rcpp_1.0.14          statmod_1.5.0        systemfonts_1.2.2    tidyselect_1.2.1    
[37] rstudioapi_0.17.1    farver_2.1.2         labeling_0.4.3       compiler_4.4.3      
[41] distributional_0.5.0

