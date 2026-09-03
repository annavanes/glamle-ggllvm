# glamle-ggllvm
Implementation code based on *Flexible latent variable models on graphs: Laplace approximated inference for multiview network data*, van Es, A., Cantoni, E., and La Vecchia, D. (2026). In preparation.

## WTO Data
The folder "data" contains the dataset used in the analysis. It consists of [WTO](https://data.wto.org) data for year 2022 and variables from `cepiigeodist` package in R, with the following columns:
- `year`: Year of the data (2022).
- `aggregation_level`: Level of aggregation for the data (3).
- `product_code`: Product code.
- `value`: Value of goods traded.
- `reporter_country`: Importer country.
- `reporter_iso3a`: Importer country code (ISO 3-letter code).
- `reporter_iso2a`: Importer country code (ISO 2-letter code).
- `partner_country`: Exporter country.
- `partner_iso3a`: Exporter country code (ISO 3-letter code).
- `partner_iso2a`: Exporter country code (ISO 2-letter code).
- `product`: Product traded.
- `reporter_region`: Importer region.
- `partner_region`: Exporter region.
- `best_orig`: Bilateral best applied simple average tariff data, see [WTO](https://data.wto.org).
- `mfn`: MFN (Most-Favoured-Nation) simple average tariff data, see [WTO](https://data.wto.org).
- `best`: Best applied if existing, replaced by mfn if missing. 
- `dist, distcap, distw, distwces, contig, comlang_off, comlang_ethno, colony, comcol, curcol, col45, smctry`: See `cepiigeodist` package description in R.

**Glossary**
Product classification: Product identification systems for classifying goods.
 
Harmonized system (HS): An international nomenclature developed by the World Customs Organization. It is arranged in six-digit codes, allowing all participating economies to classify traded goods on a common basis. Beyond the six-digit level, economies are free to introduce national distinctions for tariffs and many other purposes.
 
Multilateral Trade Negotiations (MTN): The product classification system used by the WTO for trade statistics and policy analysis. The MTN categories are defined according to the Harmonized System and consist of a two-level structure as follows:
 
22 MTN categories
72 MTN sub-categories
 
To explore MTN further, please click [here](https://stats.wto.org/Areas/TimeSeries/src/assets/WTO_Multilateral_Trade_Negotiations_Categories_2023-06-26.pdf).
 
Simple average: The unweighted average of the ad valorem or ad valorem equivalents (AVEs) of MFN applied, MFN final bound, or preference tariffs based on pre-aggregated HS six-digit averages.


## R
The folder "R" contains two folder for the Realdata analysis and the simulation codes.

### Simulations
The folder "simulation" contains the following subfolders:

#### Functions
**Functions.R**
All functions used in the simulations as well as real data analysis.

#### Poisson
**Sim_Poi_bijxij.R**

Simulation code for Model 4 with Poisson assumption , i.e. with the linear predictor $\eta_{ij}^{(k)}$ $\eta^{(k)}_{ij} =  \boldsymbol{\alpha}_{ij}^{\top} \boldsymbol{z}^{(k)} + \boldsymbol{\beta}_{ij}^{\top}\boldsymbol{x}^{(k)}_{ij}$.

**Sim_Poi_bxij.R**

Simulation code for Model 3 with Poisson assumption , i.e. with the linear predictor $\eta^{(k)}_{ij} =  \boldsymbol{\alpha}_{ij}^{\top} \boldsymbol{z}^{(k)} + \boldsymbol{\beta}^{\top}\boldsymbol{x}^{(k)}_{ij}$.

#### ZAGA
**Sim_ZAGA_bijxijwij.R**

Simulation code for Model 5 with ZAGA assumption , i.e. with the linear predictor $\eta^{(k)}_{ij} = \boldsymbol{\alpha}_{ij}^{\top} \boldsymbol{z}^{(k)} + \boldsymbol{\beta}_{ij}^{\top}\boldsymbol{x}^{(k)}_{ij} + \boldsymbol{\gamma}^{\top} \boldsymbol{w}_{ij}$.

### Realdata 
**wto_fit.R**
Fit for GLAMLE for WTO data and comparison to PPMLE fit.

## References
C. Jiang, D. La Vecchia, and R. Rastelli. GLAMLE: inference for multiview network data in the presence of latent variables, with
an application to commodities trading. Econometrics and Statistics, in press, 2025.
E. Cantoni and E. Ronchetti. Robust inference for generalized linear models. Journal of the American Statistical Association,
96(455):1022–1030, 2001.
K. Kristensen, A. Nielsen, C. Berg, H. Skaug, and B. Bell. TMB: automatic differentiation and Laplace approximation. arXiv
preprint arXiv:1509.00660, 2015.
R. A. Rigby, M. D. Stasinopoulos, G. Z. Heller, and F. De Bastiani. Distributions for modeling location, scale, and shape: Using
GAMLSS in R. Chapman and Hall/CRC, 2019.
Z. Shun and P. McCullagh. Laplace approximation of high dimensional integrals. Journal of the Royal Statistical Society: Series
B (Methodological), 57(4):749–760, 1995.