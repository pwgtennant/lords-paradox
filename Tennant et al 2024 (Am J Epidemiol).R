#### LORD'S PARADOX SIMULATION CODE FOR TENNANT ET AL 2024 ###############################
#
# Peter WG Tennant - Last updated 02 January 2024
#
# Github link: https://github.com/pwgtennant/lords-paradox
#
#### DETAILS ###############################
#
# This is the R code to accompany the research article entitled, "Lord's 'paradox' explained: The 50-year warning on the use of 'change scores' in observational data"  
# The code was prepared by Peter Tennant (University of Leeds and Alan Turing Institute) during 2022-2023 and has been updated following peer-review with the American Journal of Epidemiology in 2023-2024. 
# Peter Tennant acts as the guarantor for the accuracy of the work.  
# The code was run using R 4.3.2, but relies on several packages including CMAverse, which may need to be installed from Github.
#
# The underlying simulations are performed in Python using the DagSim package. This is called and executed within R using the reticulate package.  
# Two data structures are simulated, the first simulation (scenario1) does not include mediator-outcome confounding, the second simulation (scenario2) does include mediator-outcome confounding. 

# The simulation and analysis takes place in two stages, 
# In the first stage, two datasets are simulated for the two scenarios (with and without mediator-outcome confounding); these and used to create plots (which form Figure 3 and Supplementary Figure 1 in Tennant et al 2024).
# In the second stage, multiple datasets are simulated for the two scenarios and various models are run and their coefficients stored.
# The median coefficients and the 2.5th and 97.5th centiles for the two scenarios are then determined to form the simulation results. 
#
# NOTE: This code will create the following files in your working directory: 
# Two CSVs containing simulated datasets for the two scenarios
# A CSV containing the model coefficient estimates for all models in the two scenarios 
# Two PNG files for the two plots
#
#### REQUIRED PYTHON PACKAGES ###############################
# Please note that this code requires that python and dagsim are installed. Details of dagsim can be found at the dagsim github: https://github.com/uio-bmi/dagsim 

#### REQUIRED R PACKAGES AND FONTS ###############################

# List of R packages that are used ----
pkgs      <- c("plyr", 
               "dplyr",
               "tidyverse",
               "ggplot2",
               "gridExtra",
               "devtools",
               "showtext",
               "reticulate",
               "progress",
               "CMAverse")

# Check whether any packages are missing, and install if missing ---- 
new_pkgs  <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
if (!("CMAverse" %in% installed.packages()[,"Package"])) devtools::install_github("BS1125/CMAverse")

# Load the required packages ----  
for (itn in 1:length(pkgs)) suppressPackageStartupMessages(library(pkgs[itn], character.only = TRUE))

#Remove superfluous lists
rm(list=c("itn", "new_pkgs", "pkgs"))

# Load the font that will be used in the plots ----

font_add_google('Roboto')
showtext_auto()

#### MAKE LORD'S PARADOX FIGURES ###############################
#
# This section focuses on creating the plots for the two simulation scenarios, both in larger samples of 10000 units.

# Make sure you have the necessary python packages installed:

py_install("dagsim")
py_install("numpy")
py_install("pandas")

# Execute the simulation Python code ----

py_run_string("

import os
import sys
f = open(os.devnull, 'w')
sys.stdout = f

# Import the necessary python packages

import dagsim.base as ds
import numpy as np
import pandas as pd

#Set simulation seed
np.random.seed(1)

# SCENARIO 1 - WITHOUT MEDIATOR-OUTCOME CONFOUNDING

#Define the target path coefficients to be simulated

Path_X0_Y0  =  0.5
Path_X0_X0A =  1
Path_X0_Y1  =  0.1
Path_Y0_Y1  =  0.5
Path_X0A_X0B = 1
Path_X0B_Y1 =  0.15

#Calculate covariances (needed to work out the correct size of error terms to add to ensure variables are standardised)
Cov_X0_Y0  =  Path_X0_Y0
Cov_X0_X0B =  Path_X0_X0A*Path_X0A_X0B
Cov_Y0_X0B = (Path_X0_X0A*Path_X0A_X0B)*(Path_X0_Y0)

#Create functions to define each varaible

#1) Sex (X0)

def define_sex(Null):
  X0 = 2*np.random.binomial(n=1, p=0.5)-1
  return X0

#3) Baseline Weight (Y0)

def baseline_weight(X0):
  Y0 = (   X0*Path_X0_Y0
         + np.sqrt(1-(Path_X0_Y0**2))*np.random.normal(0, 1))
  return Y0

#4) Hall (X0A)
def hall(X0):
  X0A = (X0*1)
  return X0A

#5) Diet (X0B)
def diet(X0A):
  X0B = (X0A*1)
  return X0B

#5) Follow_up Weight (Y1)
def follow_up_weight(X0, Y0, X0B):
  Y1 =   (  X0*Path_X0_Y1
          + Y0*Path_Y0_Y1
          + X0B*Path_X0B_Y1
          + np.sqrt(1-(  Path_X0_Y1**2
                       + Path_Y0_Y1**2
                       + Path_X0B_Y1**2
                       + 2*Cov_X0_Y0*Path_X0_Y1*Path_Y0_Y1
                       + 2*Cov_X0_X0B*Path_X0_Y1*Path_X0B_Y1
                       + 2*Cov_Y0_X0B*Path_Y0_Y1*Path_X0B_Y1))*np.random.normal(0, 1))
  return Y1

#5) Weight Change (Y1_Y0)
def weight_change(Y0, Y1):
  Y1_Y0 = Y1 - Y0
  return Y1_Y0

# Create DagSim nodes for each variable based on the origin function

Null   = 0
X0     = ds.Node(name='sex',               function=define_sex,        args=[Null])
Y0     = ds.Node(name='baseline_weight',   function=baseline_weight,   args=[X0])
X0A    = ds.Node(name='hall',              function=hall,              args=[X0])
X0B    = ds.Node(name='diet',              function=diet,              args=[X0A])
Y1     = ds.Node(name='follow_up_weight',  function=follow_up_weight,  args=[X0,Y0, X0B])
Y1_Y0  = ds.Node(name='weight_change',     function=weight_change,     args=[Y0, Y1])

# Declare that these nodes belong to a DAG
DAG = ds.Graph([X0,Y0,X0A,X0B,Y1,Y1_Y0], 'DAG')

# Now simulate the data
data_dictionary = DAG.simulate(num_samples=10000, csv_name='scenario1')

# SCENARIO 2 - WITH MEDIATOR-OUTCOME CONFOUNDING

#Define the target path coefficients to be simulated

Path_X0_M0  =  0.5
Path_X0_Y0  =  0.7
Path_X0_X0A =  1
Path_X0_Y1  =  0.2
Path_M0_Y0  = -0.4
Path_M0_Y1  = -0.2
Path_Y0_Y1  =  0.5
Path_X0A_X0B = 1
Path_X0B_Y1 =  0.15

#Calculate covariances (needed to work out the correct size of error terms to add to ensure variables are standardised)
Cov_X0_M0  =  Path_X0_M0
Cov_X0_Y0  =  Path_X0_Y0+(Path_X0_M0*Path_M0_Y0)
Cov_X0_X0B =  Path_X0_X0A*Path_X0A_X0B
Cov_M0_Y0  =  Path_M0_Y0+(Path_X0_M0*Path_X0_Y0)
Cov_M0_X0B = (Path_X0_X0A*Path_X0A_X0B)*Path_X0_M0
Cov_Y0_X0B = (Path_X0_X0A*Path_X0A_X0B)*(Path_X0_Y0+(Path_X0_M0*Path_M0_Y0))


#Create functions to define each varaible

#1) Sex (X0)

def define_sex(Null):
  X0 = 2*np.random.binomial(n=1, p=0.5)-1
  return X0

#2) Physical Activity (M0)

def physical_activity(X0):
  M0 = (X0*Path_X0_M0
        + np.sqrt(1-(Path_X0_M0**2))*np.random.normal(0, 1))
  return M0

#3) Baseline Weight (Y0)

def baseline_weight(X0, M0):
  Y0 = (   X0*Path_X0_Y0
         + M0*Path_M0_Y0
         + np.sqrt(1-( Path_X0_Y0**2
                      + Path_M0_Y0**2
                      + 2*Cov_X0_M0*Path_X0_Y0*Path_M0_Y0))*np.random.normal(0, 1))
  return Y0

#4) Hall (X0A)
def hall(X0):
  X0A = (X0*1)
  return X0A

#5) Diet (X0B)
def diet(X0A):
  X0B = (X0A*1)
  return X0B

#5) Follow_up Weight (Y1)
def follow_up_weight(X0, M0, Y0, X0B):
  Y1 =   (  X0*Path_X0_Y1
          + M0*Path_M0_Y1
          + Y0*Path_Y0_Y1
          + X0B*Path_X0B_Y1
          + np.sqrt(1-(  Path_X0_Y1**2
                       + Path_M0_Y1**2
                       + Path_Y0_Y1**2
                       + Path_X0B_Y1**2
                       + 2*Cov_X0_M0*Path_X0_Y1*Path_M0_Y1
                       + 2*Cov_X0_Y0*Path_X0_Y1*Path_Y0_Y1
                       + 2*Cov_X0_X0B*Path_X0_Y1*Path_X0B_Y1
                       + 2*Cov_M0_Y0*Path_M0_Y1*Path_Y0_Y1
                       + 2*Cov_M0_X0B*Path_M0_Y1*Path_X0B_Y1
                       + 2*Cov_Y0_X0B*Path_Y0_Y1*Path_X0B_Y1))*np.random.normal(0, 1))
  return Y1

#5) Weight Change (Y1_Y0)
def weight_change(Y0, Y1):
  Y1_Y0 = Y1 - Y0
  return Y1_Y0

# Create DagSim nodes for each variable based on the origin function

Null   = 0
X0     = ds.Node(name='sex',               function=define_sex,        args=[Null])
M0     = ds.Node(name='physical_activity', function=physical_activity, args=[X0])
Y0     = ds.Node(name='baseline_weight',   function=baseline_weight,   args=[X0, M0])
X0A    = ds.Node(name='hall',              function=hall,              args=[X0])
X0B    = ds.Node(name='diet',              function=diet,              args=[X0A])
Y1     = ds.Node(name='follow_up_weight',  function=follow_up_weight,  args=[X0, M0, Y0, X0B])
Y1_Y0  = ds.Node(name='weight_change',     function=weight_change,     args=[Y0, Y1])

# Declare that these nodes belong to a DAG
DAG = ds.Graph([X0,M0,Y0,X0A,X0B,Y1,Y1_Y0], 'DAG')

# Now simulate the data
data_dictionary = DAG.simulate(num_samples=10000, csv_name='scenario2')

")

# Import the simulated data ----

scenario1              <- read.csv('scenario1.csv')
scenario2              <- read.csv('scenario2.csv')

# Transform variables to have plausible means and SDs ----
scenario1              <- scenario1 %>% 
  transform(sex = factor(sex, levels = c(-1,1), labels=c("Female", "Male"))) %>% 
  mutate(baseline_weight   =  baseline_weight*10 + 80,
         follow_up_weight  =  follow_up_weight*10 + 80,
         weight_change     =  weight_change*10)

scenario2              <- scenario2 %>% 
  transform(sex = factor(sex, levels = c(-1,1), labels=c("Female", "Male"))) %>% 
  mutate(baseline_weight   =  baseline_weight*10 + 80,
         follow_up_weight  =  follow_up_weight*10 + 80,
         physical_activity =  physical_activity*10 + 30,
         weight_change     =  weight_change*10)

# Define labels for the X (Y0) and Y (Y1) axes ----
y0_label                    <- expression(bold(paste(Y[0], ": Baseline weight (Kg)")))
y1_label                    <- expression(bold(paste(Y[1], ": Follow-up weight (Kg)")))

# The plots will be constructed from four panels, made up of two density plots of the baseline weight and follow-up weight distributions, the scatterplot of follow-up vs baseline weight, and a blank square 

# SCENARIO 1 - WITHOUT MEDIATOR-OUTCOME CONFOUNDING ----

# Draw the first density plot, for baseline weight ----

left_histogram_s1           <- ggplot(scenario1) +
                                      theme_classic(base_size = 16, 
                                                    base_family="Roboto") +
                                      theme(legend.position = "none",
                                            plot.margin = margin(-10.5, -0.2, -12, 0, "mm"),
                                            plot.background = element_rect(fill = "white"),
                                            axis.ticks = element_blank(),
                                            axis.line.x = element_blank(),
                                            axis.line.y = element_blank(),
                                            axis.text.x = element_blank(),
                                            axis.text.y = element_blank(),
                                            axis.title.x = element_blank(),
                                            axis.title.y =element_blank(),
                                            text = element_text(family = "Roboto")) +
                                      geom_density(data=subset(scenario1, sex=="Female"), 
                                                   aes(x=follow_up_weight, y=-after_stat(count)), 
                                                   bw=2.5, 
                                                   fill="#64197D", 
                                                   color="#64197D", 
                                                   linetype="dashed", 
                                                   linewidth=1, 
                                                   alpha=0.1) +
                                      geom_density(data=subset(scenario1, sex=="Male"), 
                                                   aes(x=follow_up_weight, y=-after_stat(count)),
                                                   bw=2.5, 
                                                   fill="#64197D", 
                                                   color="#64197D", 
                                                   linetype="solid", 
                                                   linewidth=1, 
                                                   alpha=0.1) +
                                      coord_cartesian(xlim=c(40,120), expand=FALSE) + 
                                      coord_flip() +
                                      geom_segment(aes(x = 40, y = 0, xend = 120, yend = 0), color="#64197D", linewidth=1) +
                                      annotate("text", x = 50, y = -50, label = "50", size = 8, colour="#64197D") +
                                      annotate("text", x = 60, y = -50, label = "60", size = 8, colour="#64197D") +
                                      annotate("text", x = 70, y = -50, label = "70", size = 8, colour="#64197D") +
                                      annotate("text", x = 80, y = -50, label = "80", size = 8, colour="#64197D") +
                                      annotate("text", x = 90, y = -50, label = "90", size = 8, colour="#64197D") +
                                      annotate("text", x = 100, y = -50, label = "100", size = 8, colour="#64197D") +
                                      annotate("text", x = 110, y = -50, label = "110", size = 8, colour="#64197D")

# Draw the second density plot, for follow-up weight ----

bottom_histogram_s1         <- ggplot(scenario1) +
                                      theme_classic(base_size = 16, 
                                                    base_family="Roboto") +
                                      theme(legend.position = "none",
                                            plot.margin = margin(-0.4, 0, 0, -3.2, "mm"),
                                            plot.background = element_rect(fill = "white"),
                                            axis.ticks = element_blank(),
                                            axis.line.x = element_blank(),
                                            axis.line.y = element_blank(),
                                            axis.text.x = element_blank(),
                                            axis.text.y = element_blank(),
                                            axis.title.x = element_blank(),
                                            axis.title.y =element_blank(),
                                            text = element_text(family = "Roboto")) +
                                      geom_density(data=subset(scenario1, sex=="Female"), 
                                                   aes(x=baseline_weight, y=-after_stat(count)), 
                                                   bounds=c(40,120), 
                                                   bw=2.5, 
                                                   fill="#C80032", 
                                                   color="#C80032", 
                                                   linetype="dashed", 
                                                   linewidth=1, 
                                                   alpha=0.1) +
                                      geom_density(data=subset(scenario1, sex=="Male"), 
                                                   aes(x=baseline_weight, y=-after_stat(count)), 
                                                   bounds=c(40,120), 
                                                   bw=2.5, 
                                                   fill="#C80032", 
                                                   color="#C80032", 
                                                   linetype="solid", 
                                                   linewidth=1, 
                                                   alpha=0.1) +
                                      coord_cartesian(xlim=c(40,120), expand=FALSE) +
                                      geom_segment(aes(x = 40, y = 0, xend = 120, yend = 0), color="#C80032", linewidth=1) +
                                      annotate("text", x = 50, y = -35, label = "50", size = 8, colour="#C80032") +
                                      annotate("text", x = 60, y = -35, label = "60", size = 8, colour="#C80032") +
                                      annotate("text", x = 70, y = -35, label = "70", size = 8, colour="#C80032") +
                                      annotate("text", x = 80, y = -35, label = "80", size = 8, colour="#C80032") +
                                      annotate("text", x = 90, y = -35, label = "90", size = 8, colour="#C80032") +
                                      annotate("text", x = 100, y = -35, label = "100", size = 8, colour="#C80032") +
                                      annotate("text", x = 110, y = -35, label = "110", size = 8, colour="#C80032")

# Draw the main scatterplot ----

scatter_s1                  <- ggplot(scenario1, aes(x = baseline_weight, y = follow_up_weight)) +
                                      theme_classic(base_size = 16, base_family="Roboto") +
                                      theme(legend.position  = "none",
                                            plot.margin      = margin(0, 0, -1.4, -3.4, "mm"),
                                            plot.background  = element_rect(fill = "white"),
                                            axis.ticks       = element_blank(),
                                            axis.line.x      = element_blank(), 
                                            axis.line.y      = element_blank(),
                                            axis.text.x      = element_blank(),
                                            axis.text.y      = element_blank(), 
                                            axis.title.x     = element_blank(),
                                            axis.title.y     = element_blank(),
                                            text             = element_text(family = "Roboto")) +
                                      geom_point(data=subset(scenario1, sex=="Female"), 
                                                 size=0.5, 
                                                 colour="#C83296", 
                                                 alpha=0.5, 
                                                 shape=16) +
                                      geom_point(data=subset(scenario1, sex=="Male"), 
                                                 size=0.5, 
                                                 colour="#009632", 
                                                 alpha=0.5, 
                                                 shape=16) +
                                      geom_abline(intercept = 0, 
                                                  slope = 1, 
                                                  color="black", 
                                                  linewidth=0.5) +
                                      stat_ellipse(data=subset(scenario1, sex=="Female"), 
                                                   geom="polygon", 
                                                   alpha=0.1, 
                                                   fill="#C83296", 
                                                   level = 0.995) +
                                      stat_ellipse(data=subset(scenario1, sex=="Male"), 
                                                   geom="polygon", 
                                                   alpha=0.1, 
                                                   fill="#009632", 
                                                   level = 0.995) +
                                      stat_ellipse(data=subset(scenario1, sex=="Female"), 
                                                   colour="#C83296", 
                                                   linetype="longdash", 
                                                   linewidth=1, 
                                                   level = 0.995) +
                                      stat_ellipse(data=subset(scenario1, sex=="Male"), 
                                                   colour="#009632", 
                                                   linetype = "solid", 
                                                   linewidth=1, 
                                                   level = 0.995) +
                                      geom_smooth(data=subset(scenario1, sex=="Female"),
                                                  method="lm", 
                                                  se = FALSE, 
                                                  colour="#C83296", 
                                                  linetype="dashed", 
                                                  fullrange=TRUE) + 
                                      geom_smooth(data=subset(scenario1, sex=="Male"),
                                                  method="lm", 
                                                  se = FALSE, 
                                                  colour="#009632", 
                                                  linetype="solid", 
                                                  fullrange=TRUE) + 
                                      geom_point(data=subset(scenario1, sex=="Female"), 
                                                 aes(x = mean(baseline_weight), y = mean(follow_up_weight)), 
                                                 size=6, 
                                                 colour="#C83296", 
                                                 shape=16) +
                                      geom_point(data=subset(scenario1, sex=="Male"), 
                                                 aes(x = mean(baseline_weight), y = mean(follow_up_weight)), 
                                                 size=6, 
                                                 colour="#009632", 
                                                 shape=16) +
                                      coord_cartesian(ylim = c(40,120), xlim=c(40,120), expand=FALSE) +
                                      annotate("text", x = 100, y = 43, label = y0_label, size = 8, fontface ="bold", colour="#C80032") + 
                                      annotate("text", x = 43, y = 100, label = y1_label, angle=90, size = 8, fontface ="bold", colour="#64197D") + 
                                      annotate("text", x = 98, y = 105, label = "Boys", size = 8, fontface ="bold", colour="#009632") + 
                                      annotate("text", x = 62, y = 56, label = "Girls", size = 8, fontface ="bold", colour="#C83296") +
                                      geom_segment(aes(x = 40, y = 40, xend = 40, yend = 120), color="#64197D", linewidth=1) +
                                      geom_segment(aes(x = 40, y = 40, xend = 120, yend = 40), color="#C80032", linewidth=1)
  
# Draw an empty plot for the final panel ----
empty                       <- ggplot()+
                                      geom_point(aes(1,1), colour="white")+
                                      theme(legend.position = "none",
                                      plot.margin = margin(0, 0, 0, 0, "mm"),
                                      axis.ticks=element_blank(), 
                                      panel.background=element_blank(), 
                                      axis.text.x=element_blank(), 
                                      axis.text.y=element_blank(),           
                                      axis.title.x=element_blank(), 
                                      axis.title.y=element_blank())

# Stitch together the final image from the four plots using grid.arrange ----

Lord_plot_s1 <- grid.arrange(left_histogram_s1, scatter_s1, empty, bottom_histogram_s1, 
                layout_matrix = matrix(c(1,2,2,2,2,2,
                                         1,2,2,2,2,2,
                                         1,2,2,2,2,2,
                                         1,2,2,2,2,2,
                                         1,2,2,2,2,2,
                                         3,4,4,4,4,4),
                                       byrow=TRUE, ncol=6))

# Save the final image as a PNG ----

ggsave(paste0("Lord_plot_scenario1", ".png"), plot = Lord_plot_s1, units = "px", height = 3000, width = 3000, dpi=300)

# SCENARIO 2 - WITH MEDIATOR-OUTCOME CONFOUNDING ----

# Draw the first density plot, for baseline weight ----

left_histogram_s2           <- ggplot(scenario2) +
                                      theme_classic(base_size = 16, 
                                                    base_family="Roboto") +
                                      theme(legend.position = "none",
                                            plot.margin = margin(-10.5, -0.2, -12, 0, "mm"),
                                            plot.background = element_rect(fill = "white"),
                                            axis.ticks = element_blank(),
                                            axis.line.x = element_blank(),
                                            axis.line.y = element_blank(),
                                            axis.text.x = element_blank(),
                                            axis.text.y = element_blank(),
                                            axis.title.x = element_blank(),
                                            axis.title.y =element_blank(),
                                            text = element_text(family = "Roboto")) +
                                      geom_density(data=subset(scenario1, sex=="Female"), 
                                                   aes(x=follow_up_weight, y=-after_stat(count)), 
                                                   bw=2.5, 
                                                   fill="#64197D", 
                                                   color="#64197D", 
                                                   linetype="dashed", 
                                                   linewidth=1, 
                                                   alpha=0.1) +
                                      geom_density(data=subset(scenario1, sex=="Male"), 
                                                   aes(x=follow_up_weight, y=-after_stat(count)),
                                                   bw=2.5, 
                                                   fill="#64197D", 
                                                   color="#64197D", 
                                                   linetype="solid", 
                                                   linewidth=1, 
                                                   alpha=0.1) +
                                      coord_cartesian(xlim=c(40,120), expand=FALSE) + 
                                      coord_flip() +
                                      geom_segment(aes(x = 40, y = 0, xend = 120, yend = 0), color="#64197D", linewidth=1) +
                                      annotate("text", x = 50, y = -50, label = "50", size = 8, colour="#64197D") +
                                      annotate("text", x = 60, y = -50, label = "60", size = 8, colour="#64197D") +
                                      annotate("text", x = 70, y = -50, label = "70", size = 8, colour="#64197D") +
                                      annotate("text", x = 80, y = -50, label = "80", size = 8, colour="#64197D") +
                                      annotate("text", x = 90, y = -50, label = "90", size = 8, colour="#64197D") +
                                      annotate("text", x = 100, y = -50, label = "100", size = 8, colour="#64197D") +
                                      annotate("text", x = 110, y = -50, label = "110", size = 8, colour="#64197D")

# Draw the second density plot, for follow-up weight ----

bottom_histogram_s2         <- ggplot(scenario2) +
                                      theme_classic(base_size = 16, 
                                                    base_family="Roboto") +
                                      theme(legend.position = "none",
                                            plot.margin = margin(-0.4, 0, 0, -3.2, "mm"),
                                            plot.background = element_rect(fill = "white"),
                                            axis.ticks = element_blank(),
                                            axis.line.x = element_blank(),
                                            axis.line.y = element_blank(),
                                            axis.text.x = element_blank(),
                                            axis.text.y = element_blank(),
                                            axis.title.x = element_blank(),
                                            axis.title.y =element_blank(),
                                            text = element_text(family = "Roboto")) +
                                      geom_density(data=subset(scenario1, sex=="Female"), 
                                                   aes(x=baseline_weight, y=-after_stat(count)), 
                                                   bounds=c(40,120), 
                                                   bw=2.5, 
                                                   fill="#C80032", 
                                                   color="#C80032", 
                                                   linetype="dashed", 
                                                   linewidth=1, 
                                                   alpha=0.1) +
                                      geom_density(data=subset(scenario1, sex=="Male"), 
                                                   aes(x=baseline_weight, y=-after_stat(count)), 
                                                   bounds=c(40,120), 
                                                   bw=2.5, 
                                                   fill="#C80032", 
                                                   color="#C80032", 
                                                   linetype="solid", 
                                                   linewidth=1, 
                                                   alpha=0.1) +
                                      coord_cartesian(xlim=c(40,120), expand=FALSE) +
                                      geom_segment(aes(x = 40, y = 0, xend = 120, yend = 0), color="#C80032", linewidth=1) +
                                      annotate("text", x = 50, y = -35, label = "50", size = 8, colour="#C80032") +
                                      annotate("text", x = 60, y = -35, label = "60", size = 8, colour="#C80032") +
                                      annotate("text", x = 70, y = -35, label = "70", size = 8, colour="#C80032") +
                                      annotate("text", x = 80, y = -35, label = "80", size = 8, colour="#C80032") +
                                      annotate("text", x = 90, y = -35, label = "90", size = 8, colour="#C80032") +
                                      annotate("text", x = 100, y = -35, label = "100", size = 8, colour="#C80032") +
                                      annotate("text", x = 110, y = -35, label = "110", size = 8, colour="#C80032")

# Draw the main scatterplot ----

scatter_s2                  <- ggplot(scenario2, aes(x = baseline_weight, y = follow_up_weight)) +
                                      theme_classic(base_size = 16, base_family="Roboto") +
                                      theme(legend.position  = "none",
                                            plot.margin      = margin(0, 0, -1.4, -3.4, "mm"),
                                            plot.background  = element_rect(fill = "white"),
                                            axis.ticks       = element_blank(),
                                            axis.line.x      = element_blank(), 
                                            axis.line.y      = element_blank(),
                                            axis.text.x      = element_blank(),
                                            axis.text.y      = element_blank(), 
                                            axis.title.x     = element_blank(),
                                            axis.title.y     = element_blank(),
                                            text             = element_text(family = "Roboto")) +
                                      geom_point(data=subset(scenario1, sex=="Female"), 
                                            size=0.5, 
                                            colour="#C83296", 
                                            alpha=0.5, 
                                            shape=16) +
                                      geom_point(data=subset(scenario1, sex=="Male"), 
                                                 size=0.5, 
                                                 colour="#009632", 
                                                 alpha=0.5, 
                                                 shape=16) +
                                      geom_abline(intercept = 0, 
                                                  slope = 1, 
                                                  color="black", 
                                                  linewidth=0.5) +
                                      stat_ellipse(data=subset(scenario1, sex=="Female"), 
                                                   geom="polygon", 
                                                   alpha=0.1, 
                                                   fill="#C83296", 
                                                   level = 0.995) +
                                      stat_ellipse(data=subset(scenario1, sex=="Male"), 
                                                   geom="polygon", 
                                                   alpha=0.1, 
                                                   fill="#009632", 
                                                   level = 0.995) +
                                      stat_ellipse(data=subset(scenario1, sex=="Female"), 
                                                   colour="#C83296", 
                                                   linetype="longdash", 
                                                   linewidth=1, 
                                                   level = 0.995) +
                                      stat_ellipse(data=subset(scenario1, sex=="Male"), 
                                                   colour="#009632", 
                                                   linetype = "solid", 
                                                   linewidth=1, 
                                                   level = 0.995) +
                                      geom_smooth(data=subset(scenario1, sex=="Female"),
                                                  method="lm", 
                                                  se = FALSE, 
                                                  colour="#C83296", 
                                                  linetype="dashed", 
                                                  fullrange=TRUE) + 
                                      geom_smooth(data=subset(scenario1, sex=="Male"),
                                                  method="lm", 
                                                  se = FALSE, 
                                                  colour="#009632", 
                                                  linetype="solid", 
                                                  fullrange=TRUE) + 
                                      geom_point(data=subset(scenario1, sex=="Female"), 
                                                 aes(x = mean(baseline_weight), y = mean(follow_up_weight)), 
                                                 size=6, 
                                                 colour="#C83296", 
                                                 shape=16) +
                                      geom_point(data=subset(scenario1, sex=="Male"), 
                                                 aes(x = mean(baseline_weight), y = mean(follow_up_weight)), 
                                                 size=6, 
                                                 colour="#009632", 
                                                 shape=16) +
                                      coord_cartesian(ylim = c(40,120), xlim=c(40,120), expand=FALSE) +
                                      annotate("text", x = 100, y = 43, label = y0_label, size = 8, fontface ="bold", colour="#C80032") + 
                                      annotate("text", x = 43, y = 100, label = y1_label, angle=90, size = 8, fontface ="bold", colour="#64197D") + 
                                      annotate("text", x = 98, y = 105, label = "Boys", size = 8, fontface ="bold", colour="#009632") + 
                                      annotate("text", x = 62, y = 56, label = "Girls", size = 8, fontface ="bold", colour="#C83296") +
                                      geom_segment(aes(x = 40, y = 40, xend = 40, yend = 120), color="#64197D", linewidth=1) +
                                      geom_segment(aes(x = 40, y = 40, xend = 120, yend = 40), color="#C80032", linewidth=1)

# Stitch together the final image from the four plots using grid.arrange ----

Lord_plot_s2 <- grid.arrange(left_histogram_s2, scatter_s2, empty, bottom_histogram_s2, 
                          layout_matrix = matrix(c(1,2,2,2,2,2,
                                                   1,2,2,2,2,2,
                                                   1,2,2,2,2,2,
                                                   1,2,2,2,2,2,
                                                   1,2,2,2,2,2,
                                                   3,4,4,4,4,4),
                                                 byrow=TRUE, ncol=6))

# Save the final image as a PNG ----

ggsave(paste0("Lord_plot_scenario2", ".png"), plot = Lord_plot_s2, units = "px", height = 3000, width = 3000, dpi=300)

#### EFFECT ESTIMATES AND COEFFICIENTS ###############################
#
# This section focuses on estimating and saving the results of various models for the two scenarios
# The simulation is repeated 10,000 times to obtain simulation intervals

# Chose the number simulations to be performed ----
Nsims <- 100

# Create a simulation progress bar ----
pb <- progress_bar$new(total = Nsims, format = ":bar :percent eta: :eta")

# Set up empty dataframes to receive the estimates from each simulation, and their summary results ----

simulation_estimates        <- data.frame(scenario1mod1=numeric(0),
                                          scenario1mod2=numeric(0),
                                          scenario1mod4=numeric(0), 
                                          scenario1mod5=numeric(0),
                                          scenario2mod1=numeric(0),
                                          scenario2mod2=numeric(0),
                                          scenario2mod3=numeric(0),
                                          scenario2mod4=numeric(0), 
                                          scenario2mod5=numeric(0))

simulation_results          <- data.frame(model=character(0), 
                                          lower=numeric(0), 
                                          point=numeric(0), 
                                          upper=numeric(0))

# Run simulation NSims times ----
for (i in 1:Nsims) {

# LOOP: Conduct Python simulation ----
py_run_string("

import os
import sys
f = open(os.devnull, 'w')
sys.stdout = f

# Import the necessary python packages

import dagsim.base as ds
import numpy as np
import pandas as pd

# SCENARIO 1 - WITHOUT MEDIATOR-OUTCOME CONFOUNDING

#Define the target path coefficients to be simulated

Path_X0_Y0  =  0.5
Path_X0_X0A =  1
Path_X0_Y1  =  0.1
Path_Y0_Y1  =  0.5
Path_X0A_X0B = 1
Path_X0B_Y1 =  0.15

#Calculate covariances (needed to work out the correct size of error terms to add to ensure variables are standardised)
Cov_X0_Y0  =  Path_X0_Y0
Cov_X0_X0B =  Path_X0_X0A*Path_X0A_X0B
Cov_Y0_X0B = (Path_X0_X0A*Path_X0A_X0B)*(Path_X0_Y0)

#Create functions to define each varaible

#1) Sex (X0)

def define_sex(Null):
  X0 = 2*np.random.binomial(n=1, p=0.5)-1
  return X0

#3) Baseline Weight (Y0)

def baseline_weight(X0):
  Y0 = (   X0*Path_X0_Y0
         + np.sqrt(1-(Path_X0_Y0**2))*np.random.normal(0, 1))
  return Y0

#4) Hall (X0A)
def hall(X0):
  X0A = (X0*1)
  return X0A

#5) Diet (X0B)
def diet(X0A):
  X0B = (X0A*1)
  return X0B

#5) Follow_up Weight (Y1)
def follow_up_weight(X0, Y0, X0B):
  Y1 =   (  X0*Path_X0_Y1
          + Y0*Path_Y0_Y1
          + X0B*Path_X0B_Y1
          + np.sqrt(1-(  Path_X0_Y1**2
                       + Path_Y0_Y1**2
                       + Path_X0B_Y1**2
                       + 2*Cov_X0_Y0*Path_X0_Y1*Path_Y0_Y1
                       + 2*Cov_X0_X0B*Path_X0_Y1*Path_X0B_Y1
                       + 2*Cov_Y0_X0B*Path_Y0_Y1*Path_X0B_Y1))*np.random.normal(0, 1))
  return Y1

#5) Weight Change (Y1_Y0)
def weight_change(Y0, Y1):
  Y1_Y0 = Y1 - Y0
  return Y1_Y0

# Create DagSim nodes for each variable based on the origin function

Null   = 0
X0     = ds.Node(name='sex',               function=define_sex,        args=[Null])
Y0     = ds.Node(name='baseline_weight',   function=baseline_weight,   args=[X0])
X0A    = ds.Node(name='hall',              function=hall,              args=[X0])
X0B    = ds.Node(name='diet',              function=diet,              args=[X0A])
Y1     = ds.Node(name='follow_up_weight',  function=follow_up_weight,  args=[X0,Y0, X0B])
Y1_Y0  = ds.Node(name='weight_change',     function=weight_change,     args=[Y0, Y1])

# Declare that these nodes belong to a DAG
DAG = ds.Graph([X0,Y0,X0A,X0B,Y1,Y1_Y0], 'DAG')

# Now simulate the data
data_dictionary = DAG.simulate(num_samples=1000, csv_name='scenario1')

# SCENARIO 2 - WITH MEDIATOR-OUTCOME CONFOUNDING

#Define the target path coefficients to be simulated

Path_X0_M0  =  0.5
Path_X0_Y0  =  0.7
Path_X0_X0A =  1
Path_X0_Y1  =  0.2
Path_M0_Y0  = -0.4
Path_M0_Y1  = -0.2
Path_Y0_Y1  =  0.5
Path_X0A_X0B = 1
Path_X0B_Y1 =  0.15

#Calculate covariances (needed to work out the correct size of error terms to add to ensure variables are standardised)
Cov_X0_M0  =  Path_X0_M0
Cov_X0_Y0  =  Path_X0_Y0+(Path_X0_M0*Path_M0_Y0)
Cov_X0_X0B =  Path_X0_X0A*Path_X0A_X0B
Cov_M0_Y0  =  Path_M0_Y0+(Path_X0_M0*Path_X0_Y0)
Cov_M0_X0B = (Path_X0_X0A*Path_X0A_X0B)*Path_X0_M0
Cov_Y0_X0B = (Path_X0_X0A*Path_X0A_X0B)*(Path_X0_Y0+(Path_X0_M0*Path_M0_Y0))


#Create functions to define each varaible

#1) Sex (X0)

def define_sex(Null):
  X0 = 2*np.random.binomial(n=1, p=0.5)-1
  return X0

#2) Physical Activity (M0)

def physical_activity(X0):
  M0 = (X0*Path_X0_M0
        + np.sqrt(1-(Path_X0_M0**2))*np.random.normal(0, 1))
  return M0

#3) Baseline Weight (Y0)

def baseline_weight(X0, M0):
  Y0 = (   X0*Path_X0_Y0
         + M0*Path_M0_Y0
         + np.sqrt(1-( Path_X0_Y0**2
                      + Path_M0_Y0**2
                      + 2*Cov_X0_M0*Path_X0_Y0*Path_M0_Y0))*np.random.normal(0, 1))
  return Y0

#4) Hall (X0A)
def hall(X0):
  X0A = (X0*1)
  return X0A

#5) Diet (X0B)
def diet(X0A):
  X0B = (X0A*1)
  return X0B

#5) Follow_up Weight (Y1)
def follow_up_weight(X0, M0, Y0, X0B):
  Y1 =   (  X0*Path_X0_Y1
          + M0*Path_M0_Y1
          + Y0*Path_Y0_Y1
          + X0B*Path_X0B_Y1
          + np.sqrt(1-(  Path_X0_Y1**2
                       + Path_M0_Y1**2
                       + Path_Y0_Y1**2
                       + Path_X0B_Y1**2
                       + 2*Cov_X0_M0*Path_X0_Y1*Path_M0_Y1
                       + 2*Cov_X0_Y0*Path_X0_Y1*Path_Y0_Y1
                       + 2*Cov_X0_X0B*Path_X0_Y1*Path_X0B_Y1
                       + 2*Cov_M0_Y0*Path_M0_Y1*Path_Y0_Y1
                       + 2*Cov_M0_X0B*Path_M0_Y1*Path_X0B_Y1
                       + 2*Cov_Y0_X0B*Path_Y0_Y1*Path_X0B_Y1))*np.random.normal(0, 1))
  return Y1

#5) Weight Change (Y1_Y0)
def weight_change(Y0, Y1):
  Y1_Y0 = Y1 - Y0
  return Y1_Y0

# Create DagSim nodes for each variable based on the origin function

Null   = 0
X0     = ds.Node(name='sex',               function=define_sex,        args=[Null])
M0     = ds.Node(name='physical_activity', function=physical_activity, args=[X0])
Y0     = ds.Node(name='baseline_weight',   function=baseline_weight,   args=[X0, M0])
X0A    = ds.Node(name='hall',              function=hall,              args=[X0])
X0B    = ds.Node(name='diet',              function=diet,              args=[X0A])
Y1     = ds.Node(name='follow_up_weight',  function=follow_up_weight,  args=[X0, M0, Y0, X0B])
Y1_Y0  = ds.Node(name='weight_change',     function=weight_change,     args=[Y0, Y1])

# Declare that these nodes belong to a DAG
DAG = ds.Graph([X0,M0,Y0,X0A,X0B,Y1,Y1_Y0], 'DAG')

# Now simulate the data
data_dictionary = DAG.simulate(num_samples=1000, csv_name='scenario2')

")

# LOOP: Import simulated data ----

scenario1              <- read.csv('scenario1.csv')
scenario2              <- read.csv('scenario2.csv')

# LOOP: Transform variables ----
scenario1              <- scenario1 %>% 
  transform(sex = factor(sex, levels = c(-1,1), labels=c("Female", "Male"))) %>% 
  mutate(baseline_weight   =  baseline_weight*10 + 80,
         follow_up_weight  =  follow_up_weight*10 + 80,
         weight_change     =  weight_change*10)

scenario2              <- scenario2 %>% 
  transform(sex = factor(sex, levels = c(-1,1), labels=c("Female", "Male"))) %>% 
  mutate(baseline_weight   =  baseline_weight*10 + 80,
         follow_up_weight  =  follow_up_weight*10 + 80,
         physical_activity =  physical_activity*10 + 30,
         weight_change     =  weight_change*10)

# LOOP: Run models and store estimates ----
  
# LOOP: Model 1 Analysis of change-score ----
s1mod1            <- lm(weight_change ~ sex, data = scenario1)   
s1mod1_estimate   <- s1mod1$coefficients[2]

s2mod1            <- lm(weight_change ~ sex, data = scenario2)   
s2mod1_estimate   <- s2mod1$coefficients[2]

# LOOP: Model 2: Analysis of follow-up adjusted for baseline ----
s1mod2            <- lm(follow_up_weight ~ sex + baseline_weight, data = scenario1)
s1mod2_estimate   <- s1mod2$coefficients[2]

s2mod2            <- lm(follow_up_weight ~ sex + baseline_weight, data = scenario2)
s2mod2_estimate   <- s2mod2$coefficients[2]

# LOOP: Model 3: g-computation ----
invisible(capture.output(
s2mod3            <- cmest(data = scenario2,
                         model="gformula",
                         exposure = "sex",
                         mediator = c("baseline_weight"),
                         mreg = list("linear"),
                         postc=c("physical_activity"),
                         postcreg=list("linear"),
                         outcome = "follow_up_weight", 
                         yreg = "linear",
                         EMint = FALSE,
                         a = "Male",
                         astar = "Female",
                         mval = list(80), nboot = 1)))

s2mod3_estimate   <- s2mod3$effect.pe[1]

# LOOP: Model 4: Analysis of change-score adjusted for baseline ----
s1mod4            <- lm(weight_change ~ sex + baseline_weight, data = scenario1)
s1mod4_estimate   <- s1mod4$coefficients[2]

s2mod4            <- lm(weight_change ~ sex + baseline_weight, data = scenario2)
s2mod4_estimate   <- s2mod4$coefficients[2]

# LOOP: Model 5: Analysis of follow-up without adjustment ----
s1mod5            <- lm(follow_up_weight ~ sex, data = scenario1)
s1mod5_estimate   <- s1mod5$coefficients[2]

s2mod5            <- lm(follow_up_weight ~ sex, data = scenario2)
s2mod5_estimate   <- s2mod5$coefficients[2]

# LOOP: Save local vectors into dataframe ----

simulation_estimates[nrow(simulation_estimates)+1,]  <- c(s1mod1_estimate, 
                                                          s1mod2_estimate, 
                                                          s1mod4_estimate, 
                                                          s1mod5_estimate,
                                                          s2mod1_estimate, 
                                                          s2mod2_estimate, 
                                                          s2mod3_estimate, 
                                                          s2mod4_estimate, 
                                                          s2mod5_estimate)

# Drop local vectors
rm(s1mod1_estimate, s1mod2_estimate, s1mod4_estimate, s1mod5_estimate, s2mod1_estimate, s2mod2_estimate, s2mod3_estimate, s2mod4_estimate, s2mod5_estimate)

# Display simulation progress
pb$tick()
}


# Extract point estimates and 95% simulation limits for each model ----
model_names                 <- list("s1mod1","s1mod2", "s1mod4","s1mod5", "s2mod1","s2mod2", "s2mod3","s2mod4","s2mod5")

for (j in 1:9) {
  
  centiles                             <- c(round(quantile(simulation_estimates[,j], 0.025), digits=2), round(quantile(simulation_estimates[,j], 0.5), digits=2), round(quantile(simulation_estimates[,j], 0.975),digits=2))
  simulation_results[nrow(simulation_results)+1,]  <- c(model_names[j], unname(as.list(centiles)))
  rm(centiles)
}

# View and save the results ----
View(simulation_results)
write.csv(simulation_results, "lords_paradox_simulation_results.csv")