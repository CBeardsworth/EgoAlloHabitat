# Repository for: *Is habitat selection in the wild shaped by individual-level cognitive biases in orientation strategy?*

Christine E. Beardsworth, Mark A. Whiteside, Philippa R. Laker, Ran Nathan, Yotam Orchan, Sivan Toledo, Jayden O. van Horik, Joah R. Madden

For any questions about the code please contact Christine at c.e.beardsworth@gmail.com

To use any data contained in this repository contact Joah at j.r.madden@exeter.ac.uk for permission.

In this repository, we have included a run-through of the R analysis [here] (https://cbeardsworth.github.io/Pheasant_OrientStrat_Habitat/) to show the outputs of the analysis without the need to run the code. For those that might want to run the code themselves, we have included three R scripts ([/R](https://github.com/CBeardsworth/NavigationHabitat/blob/master/R)) and their accompanying datasets ([/Data](https://github.com/CBeardsworth/NavigationHabitat/blob/master/Data)):

Cognition analysis and figs.R = Run the cognition analysis for the first section of the manuscript and create the figures. For this, the datasets mazeData.csv (the learning trials) and mazeRotationResults.csv (the probe trial) are required. 

iSSA analysis and bootstrapping.R = Run iSSA models and bootstrapping. This produces the datasets required for the next stage of analysis. For this code, the datasets habitat.grd (habitat information), atlas2018-strategy.csv (atlas data + id and strategy data for each bird) and FeederCoords2017_27700.csv (coordinates of feeder locations from 2017-2018) are required. The produced datasets are included in [/Data](https://github.com/CBeardsworth/NavigationHabitat/blob/master/Data) therefore to run subsequent analyses, this code does not need to be run. To develop this code we relied heavily on the code included in the supplementary material of [Signer et al. (2019)](<https://doi.org/10.1002/ece3.4823>) as well as an [online tutorial](<https://bsmity13.github.io/log_rss>) from Brian J. Smith for calculating log RSS.

Habitat analysis and Figs.R = Run the statistical models for the final section of the manuscript and create the figures. For this code, the datasets produced in the previous R script are required (habitatOrientation_coefs.csv and habitatOrientation_avail.csv). We have included [these datasets](https://github.com/CBeardsworth/NavigationHabitat/blob/master/Data) so users do not need to run the iSSA analysis and bootstrapping.R script themselves. 
