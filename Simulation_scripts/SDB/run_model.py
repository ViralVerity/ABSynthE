iteration_number_outside = 1000
iteration_count = -1

if iteration_count == -1:
    
    print("Successfully importing modules")
    
    from collections import defaultdict
    from multiprocessing.pool import ThreadPool
    
    import Tree_simulator as cts
    import file_functions 
    
    #from make_contact_dicts import *
    from make_contact_dicts_chiefdom import *
    
    from individual_class import *
    from case_class import *
    import distribution_functions
    import index_functions
    from epidemic_function import *
    from model_function import *
    
    #dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/"
    #results_path = "Looping models/Results/Fitting/"
    
   
    
    capped = False
    case_limit = 50000
    
    print("Defining parameters")
    
    
    popn_size = 7092142
    epidemic_length = 148
    cfr = 0.7
    sampling_percentage = 0.16
    SDB_start = 1
    SDB_success = 1 #proportion of burials successfully done 
    #Right now, it's binary as to when SDB starts and then how successful it is, but can gradually increase efficacy as time goes on by adding in a function for these two
    
    district_list = ["bo", 'bombali', 'bonthe', 'kailahun', 'kambia', 'kenema', 'koinadugu', 'kono', 'moyamba', 'portloko', 'pujehun', 'tonkolili', 'westernarearural', 'westernareaurban']
    
    ch_list = ['badjia', 'bagbo', 'bagbwe', 'baoma', 'bumpengawo', 'gbo', 'jaiama-bongor', 'nongobabullom', 'sittia', 'sogbini', 'yawbeko', 'biriwa', 'gbanti-kamaranka', 'gbendembungowahun', 'libeisaygahun', 'magbaimbandorhahun', 'makarigbanti', 'makenicity', 'pakimasabong', 'safrokolimba', 'sandaloko', 'sandatendaren', 'sellalimba', 'tambakka', 'benducha', 'bum', 'dema', 'imperi', 'jong', 'kongbora', 'kori', 'kowa', 'kpangakemo', 'kwamebaikrim', 'lowerbanta', 'ribbi', 'timdale', 'upperbanta', 'dea', 'jawie', 'kissikama', 'kissiteng', 'kissitongi', 'kpejebongre', 'kpejewest', 'luawa', 'malema', 'mandu', 'njaluahun', 'penguia', 'upperbambara', 'yawei', 'bramaia', 'gbinle-dixing', 'magbema', 'mambolo', 'masungbala', 'samu', 'sulima', 'tonkolimba', 'warawarabafodia', 'warawarayagala', 'dama', 'dodo', 'gaura', 'goramamende', 'kandulekpeama', 'kenemacity', 'koya_k', 'langrama', 'malegohun', 'niawa', 'nomo', 'nongowa', 'simbaru', 'smallbo', 'tunkia', 'wandor', 'dembeliasinkunia', 'diang', 'follosabadembelia', 'kasunko', 'malalmara', 'mongo', 'neya', 'nieni', 'sambaya', 'tane', 'yoni', 'fiama', 'gbane', 'gbanekandor', 'gbense', 'goramakono', 'kamara', 'koidu/new', 'lei', 'mafindor', 'nimikoro', 'nimiyama', 'sandor', 'soa', 'tankoro', 'toli', 'bagruwa', 'bumpeh', 'dasse', 'fakunya', 'kagboro', 'kaiyamba', 'kamajei', 'pangakrim', 'pejeh', 'sorogbema', 'sowa', 'bkm', 'buyaromende', 'dibia', 'kaffubullom', 'koya_pl', 'lokomasama', 'maforki', 'marampa', 'masimera', 'sandamagbolontor', 'tms', 'barri', 'gallinasperi', 'kpaka', 'kpanga-kabonde', 'makpele', 'malen', 'manosakrim', 'bocity', 'gbonkolenken', 'kafesimira', 'kakua', 'kalansogoia', 'kholifamabang', 'kolifarowalla', 'komboya', 'kunikebarina', 'kunikesanda', 'lugbu', 'niawalenga', 'selenga', 'tikonko', 'valunia', 'wunde', 'westernarearural', 'westernareaurban']
    
    distributions = distribution_functions.define_distributions() 

    print("Importing dictionaries")
    
    contact_structure = make_contact_dicts(dropbox_path)
    
    run_number = 10
    
    try:
        file_functions.make_directories(dropbox_path, results_path, run_number)
    
    except FileExistsError:
        pass

    R0_output, size_output, most_recent_tip_file, length_output = file_functions.make_summary_files(dropbox_path, results_path, run_number)
    
    if capped:
        run_out_summary = file_functions.prep_runout_summary(dropbox_path, results_path, run_number)


            
pool = ThreadPool(4)

print("Running infection model")

pool.map(run_model,(iteration_number_outside,))        
        
R0_output.close()
size_output.close()
length_output.close()
most_recent_tip_file.close()
                                  
if capped:
    run_out_summary.close()
