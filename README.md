# lords-paradox
This repository contains the analytic code to accompany the article entitled "Lord's 'paradox' explained: The 50-year warning on the use of 'change scores' in observational data". The code was prepared by Peter Tennant (University of Leeds and Alan Turing Institute) during 2022-2024 who acts as the guarantor. The final code was prepared and tested using R 4.3.2 and Python 3.11 on Windows 11. 

Two versions of the code are available, one to accompany the ArXiv preprint, and one to accompany the final article in the American Journal of Epidemiology. The preprint code is entitled "Tennant et al 2023 (ArXiv)" and the final article code is entitled "Tennant et al 2024 (Am J Epidemiol)".

Both versions simulate data in Python using the DagSim package, which is called and executed within R using the reticulate package. In the Am J Epidemiol version, this is integrated into a R Markdown file. The code relies on several packages including DagSim (for python) and CMAverse (for R), which will need to be fully installed for the code to correctly run. We have tried to include installation lines to the code, but you may encounter various problems with established installations, particularly if your R installation is not able to set up and edit your python environment.  

Both version will create CSV and PNG files in your working directory. For further details, see the annotation at the top of the relevant file.
