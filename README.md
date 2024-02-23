---
---
---

# Decision time in dynamical confidence models

This repository contains code and data used in the paper **Confidence is influenced by evidence accumulation time in dynamical decision models** (Hellmann, S., Zehetleitner, M., & Rausch, M. (preprint: <https://osf.io/preprints/psyarxiv/5ze8t/>). Data was previously published in Hellmann, et al. (2023), Shekhar & Rahnev (2021), and Pleskac & Busemeyer (2010).

## Structure:

-   dynConfiR-source package file (.tar.gz) (development version of the dynConfiR package used for analyses)
-   "ModelAnalysis": folder with code to produce simulations of models (double-increase pattern in 2DSD+, Ornstein-Uhlenbeck model simulations, and simulations and computations of optimal confidence)
-   "EmpiricalStudy": folder with
    -   "data" folder containing .RData files with the raw data (.csv's are available at the confidence database (<https://osf.io/s46pr/>)) for the first three experiments analysed
    -   Main_Script.R file which includes all analysis, from model fitting to prediction, data aggregation, visualization, and quantitative comparison. This file is structured in sections which can be folded and unfolded in RStudio
    -   several .RData files with saved results, which are loaded in Main_Script.R to avoid time-consuming computations (like model fitting), when checking reproducibility
    - a subfolder "SAT_Analysis", which contains the analysis of the fourth experiment
        - LineStudyData.csv the raw data published by Pleskac & Busemeyer (2010)
        - fitting_fcts a folder with helper functions to fit the models and simulate predictions (because the fitting routines for these experiments are not available in the dynConfiR package)
        - SAT_Main_Script.R calls the fitting functions in the fitting_fct folder and generates visualizations of the results
        - several .RDATA files including saved results 
-   sessionInfo.txt with the text-output of the sessionInfo() function
-   SessionInfo.RData including the output of sessionInfo() as R object SessionInfo (this can be used to install the exact package versions more easily)
-   InstallRequirements.R providing code to easily install the same package versions used for the analyses as stated in the SessionInfo files.
-   LICENCE a licence file for this material, which is published under the GNU GPL

## Usage:

-   Start R with package file in working directory (use your favorite CRAN mirror)

<!-- -->

```         
install.packages("dynWEV_0.0.3.tar.gz", type = "source", dependencies=TRUE,repos="http://some.cran.mirror")
```

-   If necessary, install required packages:

<!-- -->

```         
install.packages(c(""BayesFactor", "tidyverse", "Rmisc", "ggh4x", "ggpubr", "FNN")) 
```

-   To redo the whole analyses, run `EmpiricalStudy/Main_Script.R` and `EmpiricalStudy/SAT_Analysis/SAT_Main_Script.R`
-   If some `.RData` files from the repository are missing, the respective computations will be carried out again, which may take some time

### References

Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence and response time in visual perception. Psychological Review. 2023 Mar 13. <https://doi.org/10.1037/rev0000411>. Epub ahead of print. PMID: 36913292.

Shekhar, M., & Rahnev, D. (2021). The nature of metacognitive inefficiency in perceptual decision making. Psychological Review, 128(1), 45--70. <https://doi.org/10.1037/rev0000249>

Pleskac, T. J., & Busemeyer, J. (2010). Two-stage Dynamic Signal Detection: A Theory of Choice, Decision Time, and Confidence. Psychological Review, 117(3). <https://doi.org/10.1037/a0019737>

## Contact

For comments, remarks, and questions please contact me: [sebastian.hellmann\@tum.de](mailto:sebastian.hellmann@tum.de){.email}
