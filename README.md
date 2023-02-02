# lords-paradox
This is the R code to accompany the research article entitled, "Lord's 'paradox' explained: The 50-year warning on the use of 'change scores' in observational data". The code was prepared by Peter Tennant (University of Leeds and Alan Turing Institute) during 2022-2023. He acts as the guarantor for the accuracy of the work.

The code was run using R 4.2.2, but relies on several packages including CMAverse, which may need to be installed from Github. The underlying simulations are performed in Python using the DagSim package. This is called and executed within R using the reticulate package.  

The analysis broadly involves two stages, 
First a single dataset is simulated and used to create the Lord's Paradox plot.
Next, multiple datasets are simulated. Within each, various models are run and their coefficients stored.
The median coefficients and the 2.5th and 97.5th centiles are then determined to form the simulation results. 

NOTE: This code will create a CSV of a simulated dataset, a CSV containing the results, and a PNG of the Lord's Paradox plot in the working director
