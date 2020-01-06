import os

def make_directories(dropbox_path, results_path, run_number):
    
    os.mkdir(dropbox_path + results_path + str(run_number))
    os.mkdir(dropbox_path + results_path + str(run_number) + "/log_files")
    os.mkdir(dropbox_path + results_path + str(run_number) + "/trees") 
    os.mkdir(dropbox_path + results_path + str(run_number) + "/dist_mvmt")
    os.mkdir(dropbox_path + results_path + str(run_number) + "/skylines")
    os.mkdir(dropbox_path + results_path + str(run_number) + "/issues")

    os.mkdir(dropbox_path + results_path + str(run_number) + "/issues/removals")
    os.mkdir(dropbox_path + results_path + str(run_number) + "/issues/zerotaus")
    
def make_summary_files(dropbox_path, results_path, run_number):
    
    R0_output = open(dropbox_path + results_path + str(run_number) + "/R0_run.csv", 'w')
    size_output = open(dropbox_path + results_path + str(run_number) + "/epidemic_size.csv", 'w')
    most_recent_tip_file = open(dropbox_path + results_path + str(run_number) + "/most_recent_dates.csv",'w')
    length_output = open(dropbox_path + results_path + str(run_number) + "/persistence.csv",'w')
            
    most_recent_tip_file.write("number" + "," + "most_recently_sampled_tip" + "\n")
    size_output.write("number" + "," + "size" + "," + "districts_involved" + "," + "communities_involved" + "\n")
    length_output.write("number" + "," + "length_of_epidemic" + "\n")
    
    return R0_output, size_output, most_recent_tip_file, length_output


def prep_other_files(dropbox_path, results_path, run_number, iteration_count):

    tree_file = open(dropbox_path + results_path + str(run_number) + "/trees/tree_for_" + str(iteration_count) + ".txt", 'w')

    district_mvmt_file = open(dropbox_path + results_path + str(run_number) + "/dist_mvmt/mvmt_for_" + str(iteration_count) + ".csv", 'w')

    skyline_file = open(dropbox_path + results_path + str(run_number) + "/skylines/skyline_for_" + str(iteration_count) + ".csv", 'w')


    skyline_file.write("number" + ","
                       + "start_interval" + ","
                       + "end_interval" + ","
                       + "logpopize" + "\n")

    district_mvmt_file.write("DistrictOne" + ","
                             + "DistrictTwo" + ","
                             + "Times" + "\n")


    
    
    return tree_file, district_mvmt_file, skyline_file

def prep_runout_summary(dropbox_path, results_path, run_number):
    
    run_out = open(dropbox_path + results_path + str(run_number) + "/capped_epidemics.csv", 'w')
    run_out.write("Iteration number" + "," + "case_number" + "\n")
    
    return run_out

#At the start, and then also at the end if we need it if epidemic is capped
def prep_info_file(dropbox_path, results_path, run_number, index_case_individual, iteration_count):
         
    info_file = open(dropbox_path + results_path + str(run_number) + "/log_files/information_file_for_" + str(iteration_count) + ".csv", 'w')
    
    info_file.write("Individual" + "," 
                    + "Parent" + "," 
                    + "Household" + ","
                    + "District" + ","
                    + "Day_infected" + ","
                    + "Day_onset" + "," 
                    + "Day_sampled" +
                    "\n")


    info_file.write(str(index_case_individual.unique_id) + "," 
                    + "NA" + ","
                    + str(index_case_individual.hh) + "," 
                    + str(index_case_individual.dist) + ","
                    + str(0) + ","
                    + str(index_case_individual.incubation_day) + ","
                    + str(index_case_individual.incubation_day) +
                    "\n")
    
    
    return info_file



def prep_runout_file(dropbox_path, results_path, run_number, iteration_count):
    
    runout_file = open(dropbox_path + results_path + str(run_number) + "/log_files/information_file_for_" + str(iteration_count) + ".csv", 'w')
    
    runout_file.write("Individual" + "," 
        + "Parent" + "," 
        + "Household" + ","
        + "District" + ","
        + "Day_infected" + ","
        + "Day_onset" + "," 
        + "Day_sampled" +
        "\n")

    return runout_file
    