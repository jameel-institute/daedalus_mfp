# Mitigating Future Respiratory Pandemics

This repository contains the source and replication code for the paper 'Mitigating Future Respiratory Pandemics in Low-, Middle- and High-Income Countries: A Modelling Study of Health, Economic and Educational Losses', available at https://doi.org/10.21203/rs.3.rs-9418446/v1. The Daedalus epidemiological–economic model, which generates the simulation output, is written in MATLAB but the figures are plotted using Rstudio and R.

---

## 1. System Requirements

### Operating System

* Tested on Windows 11
* No non-standard hardware required

### Software Dependencies

* MATLAB R2024a (MathWorks)
* R version 4.4.1
* RStudio 2024.12.0

### MATLAB Toolboxes

* Statistics and Machine Learning Toolbox 24.1
* Parallel Computing Toolbox 24.1

### R Packages

* dplyr (1.1.4)
* purrr (1.0.2)
* tidyr (1.3.1)
* readr (2.1.5)
* stringr (1.5.1)
* fBasics (4052.98)
* fitdistrplus (1.2-2)
* forecast (8.23.0)
* scam (1.2-18)
* ggplot2 (3.5.2)
* ggh4x (0.3.0)
* ggdensity (1.0.0)
* cowplot (1.2.0)
* ggpattern (1.1.3)
* patchwork (1.3.1)

---

## 2. Installation Guide

### Instructions

1. Clone or download this repository to your local machine
2. Open MATLAB and set the working directory to the root folder of this repo
3. Open RStudio and set the working directory to the `/plotting` folder within the root

### Typical Installation Time

Less than 5 minutes on a standard desktop computer

---

## 3. Demo

### Run the Daedalus Model (MATLAB)

The paper manuscript used 5000 sample synthetic countries per income-group archetype, which takes longer to run, so in the interest of time this demo is set to 50 sample sythetic countries. 

Open MATLAB, navigate to the root directory, and run:

```
run('dd_main.m')
```

This will:

* Load input data from `/input/country_data.csv` and the disease parameter files
* Run the model 
* Save outputs to `/output/archetypes/`

### Generate Figures (RStudio)

The figures should be plotted using the model output with 5000 samples (see below), but for the demo, previously computed output is available in the files `/output/archetypes/llmic_output.zip`, `/output/archetypes/umic_output.zip` and `/output/archetypes/hic_output.zip`, which should be unzipped.

Open RStudio, navigate to the `/plotting` directory, and run:

```
source("figure_2.R")
```

et cetera. This will:

* Load the data from `/output/archetypes/`
* Generate figures corresponding to those in the manuscript
* Save figures to the same directory

### Expected Output

* Income-group archetype data underlying simulations `/output/archetypes/LLMIC_data.csv`, `/output/archetypes/UMIC_data.csv` and `/output/archetypes/HIC_data.csv`
* Income-group-, disease- and closure-strategy-specific simulation output data, e.g., `/output/archetypes/LLMIC_Influenza 2009_No Closures.csv`
* Figures: `/plotting/*.png`

### Expected Run Time

* Daedalus model: The demo with 50 sample sythetic countries took less than 1 hour to run on 8 cores.
* R plotting: The figures take up to 15m to plot.

---

## 4. Instructions for Use

To run the model on the full 5000 sample synthetic countries per income-group archetype, open MATLAB, navigate to the root directory, and change line 14 to:

```
nsamples = 5000;
```

Then run:

```
run('dd_main.m')
```

as for the demo (above). 

### Expected Run Time

* Daedalus model: The full 5000 sample synthetic countries per income-group archetype took approx. 16 hours on 32 cores.

---

## License

This project is licensed under the MIT License.

---

## Code Availability

The code is publicly available in this repository.

---

## Notes

* All file paths are relative and assume the repository structure described above