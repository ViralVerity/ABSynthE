import os
import sys
import yaml
import datetime as dt


def parse_population_information(config_file):

    with open(config_file,"r") as f:
        try:
            population_config = yaml.load(f, Loader=yaml.FullLoader)
        except:
            sys.stderr.write(f'Error: failed to read config file. Ensure your file in correct yaml format.\n')
            sys.exit(-1)
    
    if type(population_config) != dict:
        sys.stderr.write("Config file is not correctly formatted.\n")
        sys.exit(-1)
        
    return population_config


def make_directories(config):

    if not config["output_directory"]:
        config["output_directory"] = f'absynthe_results_{dt.date.today()}'
    
    if not os.path.exists(config["output_directory"]):
        os.mkdir(config["output_directory"])

    if not os.path.exists(os.path.join(config["output_directory"], "log_files")):
        os.mkdir(os.path.join(config["output_directory"], "log_files"))
    elif not config["overwrite"]:
        for f in  os.listdir(os.path.join(config["output_directory"], "log_files")):
            if not f.startswith("."):
                sys.stderr.write(f"Results already exist at {config['output_directory']}. Please provide a different output directory or use the --overwrite flag.\n")
                sys.exit(-1)

    if not os.path.exists(os.path.join(config["output_directory"],"trees")):
        os.mkdir(os.path.join(config["output_directory"],"trees"))
        os.mkdir(os.path.join(config["output_directory"],"dist_mvmt"))
        os.mkdir(os.path.join(config["output_directory"],"ch_mvmt"))
        os.mkdir(os.path.join(config["output_directory"], "skylines"))
        os.mkdir(os.path.join(config["output_directory"], 'lineages'))

    return config
        
def make_summary_files(config):
    
    if config["calculate_R0"]:
        R0_output = open(os.path.join(config["output_directory"],"R0_run.csv"), 'w')
    else:
        R0_output = None

    config["R0_output"] = R0_output
        
    config["size_output"] = open(os.path.join(config["output_directory"], "epidemic_size.csv"),'w')
    config["length_output"] = open(os.path.join(config["output_directory"], "persistence.csv"), 'w')
    
    if config["output_tree"] or config["calculate_R0"] or config["output_ltt"] or config["output_skyline"]:
        config["most_recent_tip_file"] = open(os.path.join(config["output_directory"], "most_recent_dates.csv"),'w')
        config["most_recent_tip_file"].write("number,most_recently_sampled_tip\n")
    
    config["size_output"].write("number,size,districts_involved,chiefdoms_involved\n")
    config["length_output"].write("number,length_of_epidemic\n")
    
    return config

def prep_movement_files(output_directory, iteration_count):

    district_mvmt_file = open(os.path.join(output_directory,"dist_mvmt", f'mvmt_for_{iteration_count}.csv'), 'w')
    ch_mvmt_file = open(os.path.join(output_directory,"ch_mvmt",f"mvmt_for_{iteration_count}.csv"), 'w')
    
    district_mvmt_file.write("district_one,district_two,times\n") 
    ch_mvmt_file.write("chiefdom_one,chiefdom_two,times\n")                      

    return district_mvmt_file, ch_mvmt_file

def prep_runout_summary(output_directory):
    
    run_out = open(os.path.join(output_directory, "capped_epidemics.csv"), 'w')
    run_out.write("Iteration number,case_number\n")
    
    return run_out

#At the start, and then also at the end if we need it if epidemic is capped
def prep_info_file(output_directory, index_case_individual, iteration_count):
         
    info_file = open(os.path.join(output_directory,"log_files",f"information_file_for_{iteration_count}.csv"), 'w')
    info_file.write("Individual,Parent,Household,Chiefdom,District,Day_infected,Day_onset,Day_sampled\n") 
    info_file.write(f"{index_case_individual.unique_id},NA,{index_case_individual.hh},{index_case_individual.ch},{index_case_individual.dist},0,{index_case_individual.incubation_day},{index_case_individual.incubation_day}\n")

    return info_file

    