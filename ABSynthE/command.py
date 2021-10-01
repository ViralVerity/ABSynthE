#!/usr/bin/env python3
from absynthe import __version__

import argparse

from collections import defaultdict
from multiprocessing.pool import ThreadPool

import absynthe.stochastic.tree_simulator as tree_sim #was called cts before
from absynthe.stochastic.epidemic_functions import *

import absynthe.set_up.file_functions as file_funcs
from absynthe.set_up.make_contact_dicts import *
import absynthe.set_up.index_functions as index_funcs
import absynthe.set_up.distribution_functions as dist_funcs

from abysnthe.classes.individual_class import *
from abysnthe.classes.case_class import *



def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False,
    description=misc.preamble(__version__))

    parser.add_argument("--input-directory", "-indir", dest="input_directory", help="directory containing contact structure jsons")
    parser.add_argument("--output-directory", "-outdir", dest="output_directory")
    parser.add_argument("--case-limit", dest="case_limit", help="Cap the epidemic at this many cases")
    parser.add_argument("--day-limit", dest="day_limit", help="Cap the epidemic at this many days.")
    parser.add_argument("--cfr", help="Set the case fatality rate as a number between 0 and 1. Default is 0.7 for Ebola", default=0.7)
    parser.add_argument("--sampling-percentage", "-spct",  dest="sampling_percentage", help="Percentage of cases sampled, used in generating the phylogeny from cases. Default is 0.16 for Ebola", default=0.16)

    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)
    
    input_directory = args.input_directory
    output_directory = args.output_directory
    case_limit = args.case_limit
    day_limit = args.day_limit
    cfr = args.cfr
    args.sampling_percentage = args.sampling_percentage
    
    cwd = os.getcwd()
    thisdir = os.path.abspath(os.path.dirname(__file__))
    
    epidemic_length = 148 #can't remember if this is used when capped is false, so leaving in here for a minute
    
    distributions = dist_funcs.define_distributions() #this is ebola specific, so would be nice here to have user input
    contact_structure = make_contact_dicts(input_directory)
                
    file_funcs.make_directories(output_directory)
    R0_output, size_output, most_recent_tip_file, length_output = file_funcs.make_summary_files(output_directory)
    
    if case_limit or day_limit:
        run_out_summary = file_functions.prep_runout_summary(output_directory) #need to get those secondary args
    
    #where does the info file get prepped?


    pool = ThreadPool(25)

    pool.map(run_model,(iteration_number_outside,)) #calls the run_model.py script       
            
    R0_output.close()
    size_output.close()
    length_output.close()
    most_recent_tip_file.close()
                                    
    if case_limit or day_limit:
        run_out_summary.close()
