# ABSynth

Agent-based simulator for Ebola Virus

## Simulation scripts

Contains scripts for running various different versions of the simulator, and associated modules.

"Simulate_epidemic" - Runs the simulator, with options for number of iterations; capping/no capping etc. Imports all of the following modules.
"case_class.py" - Defines the case class and associated methods
"distribution_functions.py" - defines cdfs for infection parameters eg likelihood of death over time
"epidemic_function.py" - runs one iteration of the epidemic through the defined contact structure
"file_functions.py" - Preps results files
"index_functions.py" - Makes index case for epidemic
"individual_class.py" - Defines individual class and associated methods
"make_contact_dicts.py" - makes dictionaries of various aspects of the gross contact structure 
NB the data is currently not on this repo as the files are too large to be uploaded



Non-stoch model - to generate the same trees each time for testing other hypotheses

Tree_simulator - generates the coalescent tree and corresponding skyline stochastically 
Tree simulator nonstoch - the same but non-stochastically


## Results scripts

Notebooks for analysing results and producing figures from simulation runs

## To do:

- Tidy up the results scripts - what does each notebook do, what code is useful in each one

- Optimise the tree_simulator: lprof and stuff
(can probably get rid of some of the loops maybe?)

- The infection simulator could do with tidying up of which variables are actually needed and used - lots of dictionaries and things around  ###DONE
- Get rid of try/except blocks if possible - replace with if statements
- Go through all the if else statement in the main loop to make sure that they're not duplicating ###DONE

- Make node one class with flags to say transmission/normal/coalescent
- Same for tree

- Remove day cap from the normal model, add the adding days bit #### DONE
- Add an if statement to check whether I want to run the whole epidemic or just the first section, then can add it to the runtime exception to break out if we don’t want the whole epidemic #### DONE

- Check that every day is calculated in the day dict - which bits are inclusive and exclusive

- Assign functions to classes ### DONE


New functions in main simulator:
- Make files (eg with headers etc) #### DONE
- Write to file #### Not actually needed because short now
- Organise dictionaries (ie import and make agent location dict) Not sure. Might make it easier to read #### DONE
- Index case initialising ### DONE
- Remove cases ###Don't think it's really necessary
- Get output - write tree file etc, get the tree ### Don't think it's really necessary





