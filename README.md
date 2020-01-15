# ABSynthE

Agent-based simulator for Ebola Virus.

Still very much a work in progress

## Simulation scripts

Contains scripts for running various different versions of the simulator, and associated modules.

"Simulate_epidemic" - Runs the simulator, with options for number of iterations; capping/no capping etc. Imports all of the following modules:
- "case_class.py" - Defines the case class and associated methods
- "distribution_functions.py" - defines cdfs for infection parameters eg likelihood of death over time
- "epidemic_function.py" - runs one iteration of the epidemic through the defined contact structure
- "file_functions.py" - Preps results files
- "index_functions.py" - Makes index case for epidemic
- "individual_class.py" - Defines individual class and associated methods
- "make_contact_dicts.py" - makes dictionaries of various aspects of the gross contact structure 

NB the contact data is currently not on this repo as the files are too large to be uploaded. It is therefore tricky to run the simulator



- "Tree_simulator" - generates the coalescent tree and corresponding skyline stochastically. Contains some functions, and also runs:
- "tree_class.py": - contains tree class definition with associated methods
- "node_class.py" - contains node class definition and associated methods


Non-stoch model - to generate the same trees each time for testing other hypotheses

Tree simulator nonstoch - the same but non-stochastically

Testing_notebook - fairly messy jupyter for testing individual functions and exploring package functionality


## Results scripts

Notebooks for analysing results and producing figures from simulation runs







