## Run analysis, write model results

## Before: empty
## After: folder with the assessment run and the retros



# Authors: Rosana Ourens (adapted from Laura Wise for pil.27.8c9a) 

# Date: 2026

# Load libraries ----------------------------------------------------------
library(r4ss)
library(icesTAF)

# Working directory and folders -------------------------------------------

# check working directory

getwd()


# Run script for stock assessment -----------------------------------------

#*************************************************
#*      RUN THE MODEL in SS3
#*************************************************

#directories 
wd<-getwd()
run_dir <- paste0(wd,"/model/run")
SS3_exe_path <- paste0(wd,"/boot/software/ss3_win.exe")
# create output folder and sub-folders using the function mkdir from icesTAF
mkdir("model/run")


# Get input files --------------------------------------------------------

cp("boot/data/data.ss", "model/run")
cp("boot/data/control.ss", "model/run")
cp("boot/data/forecast.ss", "model/run")
cp("boot/data/starter.ss", "model/run")
cp("boot/data/wtatage.ss", "model/run")


# Run SS from r 

r4ss::run(dir = run_dir, exe = SS3_exe_path,show_in_console = T,skipfinished=F)



# Session info ------------------------------------------------------------

sessionInfo()

# End of script -----------------------------------------------------------

rm(list=ls())


