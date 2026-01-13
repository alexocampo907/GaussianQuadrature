# Calculating True Values in Causal Inference Simulation Studies via Gaussian Quadrature

This repository contains the reproduction code for the manuscript:  
**"Revealing the truth: Calculating True Values in Causal Inference Simulation Studies via Gaussian Quadrature"**

---

## 📋 General Information

- **Software:** R version 4.4.3
- **Data Requirements:** None. All data are simulated within the scripts. 
  - *Note:* Simulation results can alternatively be loaded from `Results/R Objects/`.
- **Hardware Note:** Computation times estimated below are based on a standard workstation.

### Configuration
No manual edits are required to run the primary analysis. However, parameters can be modified at the top of the scripts:
- `load.sim = F`: Change to `T` to load pre-calculated intermediate results (`.rds` files) instead of re-running simulations.
- `iter <- 10000`: Reduce this number to decrease runtime for testing purposes.

---

## 📂 File Structure & Contents

The R scripts below reproduce the specific sections and figures cited in the paper. Each `.R` file can be run in isolation.

### 1. `3_Confounding.R`
* **Description:** Reproduces the results and figures presented in Section 3 regarding confounding bias and integration.
* **Output:** Reproduces **Figures 4-8**.
* **Estimated Runtime:** ~100 minutes.
    * *Tip:* To reduce runtime, lower the iterations variable (`iter <- 10000`).

### 2. `4.1_CDE_continuous.R`
* **Description:** Performs Monte Carlo (MC) integration and Gauss-Hermite quadrature for the Controlled Direct Effect (CDE) with a continuous endpoint.
* **Output:** Reproduces the **top panel of Figure 10**.
* **Estimated Runtime:** ~11 minutes.
    * *Tip:* To reduce runtime, lower the iterations variable (`iter <- 10000`).

### 3. `4.1_CDE_binary.R`
* **Description:** Performs MC integration and Gauss-Hermite quadrature for the CDE with a binary endpoint.
* **Output:** Reproduces the **bottom panel of Figure 10**.
* **Estimated Runtime:** ~13 minutes.
    * *Tip:* To reduce runtime, lower the iterations variable (`iter <- 10000`).

### 4. `4.2_RMST.R`
* **Description:** Performs MC integration and Gauss-Hermite quadrature for mediation effects on the Restricted Mean Survival Time (RMST) scale.
* **Output:** Reproduces **all panels of Figure 11**.
* **Estimated Runtime:** ~8 minutes.
    * *Tip:* To reduce runtime, lower the iterations variable (`iter <- 10000`).

### 5. `Results/` Directory
A local directory where all figures and numerical outputs are automatically saved upon execution of the scripts.
* `Results/Plots/`: Stores generated figures.
* `Results/R Objects/`: Stores simulation results (intermediate `.rds` files).

---

## 💻 Computational Environment & Dependencies

The primary analysis was conducted using **R version 4.4.3** on Ubuntu 22.04.5 LTS.

**Key Methodological Packages:**
* `mvQuad`
* `spatstat.random` (for quadrature rules)

<details>
<summary><strong>Click to expand full Session Info</strong></summary>

```r
R version 4.4.3 (2025-02-28)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

Attached base packages:
[1] stats      graphics   grDevices  utils      datasets   methods    base     

Other attached packages:
 [1] gridExtra_2.3        spatstat.random_3.3-3 spatstat.geom_3.3-6   spatstat.univar_3
