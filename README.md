# ABSynth

Agent-based simulator for Ebola Virus

## Simulation scripts

Contains scripts for running various different versions of the simulator.

Computer model and server model are functionally the same, but have different file pathways depending on whether they are being run locally or on the server.

No caps model - no caps on days or case numbers

#### Above are to be consolidated, so that a command line option will tell it whether to run with caps or not

Non-stoch model - to generate the same trees each time for testing other hypotheses

Tree_simulator - generates the coalescent tree and corresponding skyline stochastically 
Tree simulator nonstoch - the same but non-stochastically

## Results scripts

Notebooks for analysing results and producing figures from simulation runs

## To do:

- Tidy up the results scripts - what does each notebook do, what code is useful in each one

- Optimise the tree_simulator: lprof and stuff
(can probably get rid of some of the loops maybe?)

- The infection simulator could do with tidying up of which variables are actually needed and used - lots of dictionaries and things around
- Get rid of try/except blocks if possible - replace with if statements
- Go through all the if else statement in the main loop to make sure that they're not duplicating
- Make node one class with flags to say transmission/normal/coalescent
- Same for tree

- Remove day cap from the normal model, add the adding days bit
- Add an if statement to check whether I want to run the whole epidemic or just the first section, then can add it to the runtime exception to break out if we don’t want the whole epidemic

- Check that every day is calculated in the day dict - which bits are inclusive and exclusive


New functions in main simulator:
- Make files (eg with headers etc)
- Write to file
- Organise dictionaries (ie import and make agent location dict) Not sure. Might make it easier to read
- Index case initialising
- Remove cases
- Get output - write tree file etc, get the tree





