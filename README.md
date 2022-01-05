# ABSynthE

**A**gent **B**ased **Synth**etic **E**pidemic

<img src="./logo/ABSynthE_logo.png" width="350">


## absynthe

Installable tool to run an ebola epidemic through the population of Sierra Leone in 2014. 
It can also simulate a coalescent tree with accompanying epidemic statistics.

Detailed information about the model algorithm can be found in "ABSynthE description.docx" 

To install, run: 

`python setup.py install`

To set up a general run:

`absynthe --input-directory PATH/TO/CONTACT/STRUCTURE --population-config PATH/TO/POPULATION/CONFIG`

Which will run with the defaults (currently set to Ebola in SLE defaults eg sampling percentage is 0.16, case fatality is 0.7).
This will output log files for the epidemic (i.e. information about every case, length, number of cases).

To simulate a coalescent tree, use the flag `--output-tree`

To see all options available run:

`absynthe --help`


## SLE_EBOV_input_files

Contains population structure of Sierra Leone in 2014 based on the census that came out that year. Also contains a config file with lists of the names of different contact levels eg districts and chiefdoms

To run ABSynthE on a different population, these files will all need to be present.

## Fitting 

Contains observed data from the 2013-2016 Ebola epidemic in Sierra Leone to fit the model and obtain values for transmission parameters 

Also contains scripts for running fitting using pyabc (https://pyabc.readthedocs.io/en/latest/index.html) and therefore calculating summary stats for simulated data

## Results scripts

Notebooks for analysing results and producing figures from simulation runs

- "Analyse single file" gives spread, persistence and some network analysis of individual runs
- "Comparing skygrid and skyline" draws a skygrid from a BEAST log file, plots a simulated skyline on the same axis, and compares using regression and KDE.
- "Getting tree from log file" for when a tree was not simulated but a complete log file is available
- "Summary results" gives summary results from many runs of the simulator eg distribution of epidemic lengths





