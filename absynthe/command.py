#!/usr/bin/env python3
from absynthe import __version__

import argparse

from collections import defaultdict
from multiprocessing.pool import ThreadPool

from absynthe.stochastic.run_model import *

import absynthe.set_up.file_functions as file_funcs
from absynthe.set_up.make_contact_dicts import *
import absynthe.set_up.index_functions as index_funcs
import absynthe.set_up.distribution_functions as dist_funcs

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False,
    description=preamble(__version__))

    parser.add_argument("--input-directory", "-indir", dest="input_directory", help="directory containing contact structure jsons")
    parser.add_argument("--output-directory", "-outdir", dest="output_directory")
    parser.add_argument("--population-config", "-popcon", dest="population_config")
    parser.add_argument("--case-limit", dest="case_limit",type=int, help="Cap the epidemic at this many cases")
    parser.add_argument("--day-limit", dest="day_limit", type=int,help="Cap the epidemic at this many days.")
    parser.add_argument("--cfr", type=float, help="Set the case fatality rate as a number between 0 and 1. Default is 0.7 for Ebola", default=0.7)
    
    #sampling delay after symptoms can be added in - probably wants a distribution, so maybe this would just be a flag
    parser.add_argument("--sampling-percentage", "-spct",type=float, dest="sampling_percentage", help="Percentage of cases sampled, used in generating the phylogeny from cases. Default is 0.16 for Ebola", default=0.16)
    parser.add_argument("--sampling-scheme", "-scheme", dest="sampling_scheme", help="Sampling scheme to make tree. Currently only uniform over time and space is available.")

    parser.add_argument("--output-tree", dest="output_tree", action="store_true", help="Output newick string of tree when epidemic is logged")
    parser.add_argument("--output-skyline", dest="output_skyline", action="store_true",help="Output skyline when tree is generated.")
    parser.add_argument("--output-ltt", dest="output_ltt", action="store_true", help="Calculate lineages through time when tree is generated.")
    parser.add_argument("--calculate-R0", "-R0", dest="calculate_R0", action="store_true", help="Calculate R0 when outbreak is logged")
    
    parser.add_argument("--number-model-iterations",type=int, dest="number_model_iterations", help="Number of times the stochastic epidemic is run. Default is 100", default=1000)
    parser.add_argument("--log-every", dest="log_every",type=int, help="Frequency of logging epidemics in model states. Default is 10%% of number_model_iteration ", default=0.1)
    
    parser.add_argument("--overwrite", action="store_true", help="overwrite results in output directory")
    parser.add_argument("--verbose", action="store_true", help="prints more information as it's runnign")
    parser.add_argument("-h","--help",action="store_true",dest="help")
    
    parser.add_argument("--fitting", action="store_true")

    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)

    config = {}
    config["number_model_iterations"] = args.number_model_iterations
    config["log_every"] =  args.log_every*args.number_model_iterations
    config["cfr"] = args.cfr
    
    config["sampling_percentage"] = args.sampling_percentage
    config["sampling_scheme"] = "uniform"

    config["input_directory"] = args.input_directory
    config["output_directory"] = args.output_directory
    config["case_limit"] = args.case_limit
    config["day_limit"] = args.day_limit #default is 148 for ebola in SLE for fitting (ie the exponential start) - not this any more, need to change

    config["output_tree"] = args.output_tree
    config["output_skyline"] = args.output_skyline
    config["output_ltt"] = args.output_ltt
    config["calculate_R0"] = args.calculate_R0

    config["overwrite"] = args.overwrite
    config["verbose"] = args.verbose
    
    #from ABC-SMC fitting - will hard code in here, it's for ebov SLE exponential process
    config["a"] = #??
    config["b"] = #??
    config["c"] = #??

    sys.stdout.write("Setting up for running epidemics\n")
    
    cwd = os.getcwd()
    thisdir = os.path.abspath(os.path.dirname(__file__))
    
    config = file_funcs.make_directories(config)
    config = file_funcs.make_summary_files(config)

    config["distributions"] = dist_funcs.define_distributions() #this is ebola specific, so would be nice here to have user input
    config["population_structure"] = file_funcs.parse_population_information(args.population_config)
    config = make_contact_dicts(config["input_directory"], config)
    
    if config["case_limit"] or config["day_limit"]:
        config["capped"] = True
        run_out_summary = file_funcs.prep_runout_summary(config["output_directory"])
    else:
        config["capped"] = False
        run_out_summary = ""
    
    print("\n**** CONFIG ****")
    no_print = ["population_structure", "distributions"]
    for k in sorted(config):
        if k not in no_print:
            print(f" - {k}: {config[k]}")
    
    #see if the multi-threading still works
    pool = ThreadPool(25) #might need to use star map to use lots of arguments? or put in a config dict
    pool.map(run_model,(config, )) #calls the run_model.py script       
            
    config["files"]["R0_output"].close()
    config["files"]["size_output"].close()
    config["files"]["length_output"].close()
    config["files"]["most_recent_tip_file"].close()
                                    
    if case_limit or day_limit:
        run_out_summary.close()


if __name__ == '__main__':
    main()


def preamble(v):
    print(f"""\n
    
                ABSynthE: Agent Based Synthetic Epidemic
                                {v}
                            Verity Hill
                           Andrew Rambaut
                        Edinburgh University

                \n""")
