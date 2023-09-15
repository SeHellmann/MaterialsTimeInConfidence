---
---
---

# Decision time in dynamical confidence models

This repository contains code and data used in the paper **Confidence is influenced by decision time in dynamical decision models** (Hellmann, S., Zehetleitner, M., & Rausch, M. (preprint). 
Data was previously published in Hellmann, et al. (2023) and Shekhar & Rahnev (2021). 
The most recent author manuscript of the article is available here: ..... 

## Structure:

-   dynConfiR-source package file (.tar.gz) (development version of the dynConfiR package used for analyses)
-   "ModelAnalysis": folder with code to produce figures in the formal model analysis section including simulation of KL distances and mean confidence in the dynaViTE model and a shiny application to play with parameters in figure X.
-   "EmpiricalStudy": folder with
  - "data" folder containing .RData files with the raw data (.csv's are available at the confidence database (https://osf.io/s46pr/))
  - Main_Script.R file which includes all analysis, from model fitting to prediction, data aggregation, visualization, and quantitative comparison. This file is structured in sections which can be folded and unfolded in RStudio
  - several .RData files with saved results, which are loaded in Main_Script.R to avoid time-consuming computations (like model fitting), when check reproducibility
- sessionInfo.txt with the text-output of the sessionInfo() function
- SessionInfo.RData including the output of sessionInfo() as R object SessionInfo (this can be used to install the exact package versions more easily)

## Usage:

-   Start R with package file in working directory (use your favorite CRAN mirror)

<!-- -->

    install.packages("dynWEV_0.0.3.tar.gz", type = "source", dependencies=TRUE,repos="http://some.cran.mirror")

-   If necessary, install required packages:

<!-- -->

    install.packages(c(""BayesFactor", "tidyverse", "Rmisc", "ggh4x", "ggpubr", "FNN")) 
    
-   To redo the whole analyses, run `EmpiricalStudy/Main_Script.R`
-   If some `.RData` files from the repository are missing, the respective computations will be carried out again, which may take some time

### References

Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence and response time in visual perception. Psychological Review. 2023 Mar 13. <https://doi.org/10.1037/rev0000411>. Epub ahead of print. PMID: 36913292.


## Contact

For comments, remarks, and questions please contact me: [sebastian.hellmann\@ku.de](mailto:sebastian.hellmann@ku.de){.email}
