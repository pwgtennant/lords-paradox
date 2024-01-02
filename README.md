# lords-paradox
This repository contains the analytic code to accompany the article entitled "Lord's 'paradox' explained: The 50-year warning on the use of 'change scores' in observational data". 

Two versions of the code are available, one to accompany the ArXiv preprint, and one to accompany the final article in the American Journal of Epidemiology. The preprint code is entitled "Tennant et al 2023 (ArXiv)" and the final article code is entitled "Tennant et al 2024 (Am J Epidemiol)".

The code was prepared by Peter Tennant (University of Leeds and Alan Turing Institute) during 2022-2024 who acts as the guarantor. The final code was prepared and tested using R 4.3.2 on Windows 11. The code  relies on several packages including DAGSIM (for python) and CMAverse (for R), which will need to be fully installed for the code to correctly run. You may encounter problems if your R installation is not able to set up a virtual python environment and if the virtual environment does not have all the necessary packages installed.  

Both versions simulate data in Python using the DagSim package, which is called and executed within R using the reticulate package.

Both version will create CSV and PNG files in your working directory. For further details, see the annotation at the top of the relevant file.
