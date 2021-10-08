import os
import sys
import yaml


def parse_population_information(configfile):

    with open(configfile,"r") as f:
        try:
            population_config = yaml.load(f, Loader=yaml.FullLoader)
        except:
            sys.stderr.write(cyan(f'Error: failed to read config file. Ensure your file in correct yaml format.\n'))
            sys.exit(-1)
        
    return population_config


def make_directories(output_directory):
    
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    if not os.path.exists(os.path.join(output_directory, "log_files")):
        os.mkdir(os.path.join(output_directory, "log_files"))
    else:
        sys.stderr.write(f"Results already exist at {output_directory}. Please provide a different output directory.\n")
        sys.exit(-1)

    os.mkdir(os.path.join(output_directory,"trees"))
    os.mkdir(os.path.join(output_directory,"dist_mvmt"))
    os.mkdir(os.path.join(output_directory,"ch_mvmt"))
    os.mkdir(os.path.join(output_directory, "skylines"))
    os.mkdir(os.path.join(output_directory, 'lineages'))
    
def make_summary_files(output_directory):
    
    if config["calculate_R0"]:
        R0_output = open(os.path.join(output_directory,"R0_run.csv"), 'w')
    else:
        R0_output = None
        
    size_output = open(os.path.join(output_directory, "epidemic_size.csv"),'w')
    most_recent_tip_file = open(os.path.join(output_directory, "most_recent_dates.csv"),'w')
    length_output = open(os.path.join(output_directory, "persistence.csv"), 'w')
    
    most_recent_tip_file.write("number,most_recently_sampled_tip\n")
    size_output.write("number,size,districts_involved,communities_involved\n")
    length_output.write("number,length_of_epidemic\n")
    
    return R0_output, size_output, most_recent_tip_file, length_output


def prep_movement_files(output_directory, iteration_count):

    district_mvmt_file = open(os.path.join(output_directory,"dist_mvmt", f'mvmt_for_{iteration_count}.csv', 'w'))
    ch_mvmt_file = open(os.path.join(output_directory,"ch_mvmt",f"mvmt_for_{iteration_count}.csv", 'w'))
    
    district_mvmt_file.write("district_one,district_two,times\n") 
    ch_mvmt_file.write("chiefdom_one,chiefdom_two,times\n")                      

    return district_mvmt_file, ch_mvmt_file

def prep_runout_summary(output_directory):
    
    run_out = open(os.path.join(output_directory, "capped_epidemics.csv"), 'w')
    run_out.write("Iteration number,case_number\n")
    
    return run_out

#At the start, and then also at the end if we need it if epidemic is capped
def prep_info_file(output_directory, index_case_individual, iteration_count):
         
    info_file = open(os.path.join(output_directory,"log_files",f"information_file_for_{iteration_count}.csv", 'w'))
    info_file.write("Individual,Parent,Household,Chiefdom,District,Day_infected,Day_onset,Day_sampled\n") 
    info_file.write(f"{index_case_individual.unique_id},NA,{index_case_individual.hh},{index_case_individual.comm},{index_case_individual.dist},0,{index_case_individual.incubation_day},{index_case_individual.incubation_day}\n")

    return info_file

    