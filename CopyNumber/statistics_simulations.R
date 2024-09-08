# files  = (00,01,02,03)

dir_names = paste0("../",1:12)
names = paste0(1:12)

# list of average score (among the 20 repetitions of the same simulation type) 
# of the model choosen as best with model selection
my_list <- list(simulation = names,
                    MAE = c(),
                    RI = c()
                    )


original_dir=getwd()
# for loop on simulations
for (i in dir_names) {

    setwd(i)
    print(i)
    simulation_directory=getwd()
    print(simulation_directory)
    
    # list where to store the values of the best model MAe of a simulation 
    # to then perform the average of these values to put in the main list MAE
    MAE_list = c ()
    RI_list = c ()


    simulation_directories = list.dirs(path = simulation_directory, recursive = FALSE, full.names = TRUE) 
    directories_number = length(simulation_directories)
    print(paste0("number of directories for each simulation",directories_number))

    # for loop on iteration of simulation
    for (i in 1:directories_number){ 

        setwd(paste0("simulation_iteration_",i,"/results"))
        single_directory = getwd()
        # sapply(BiblioDir,function(dir){length(list.files(dir,pattern='pdf'))})

        all_sim = readRDS("all_sim_input_prepare_input_data.rds")
        S = length(unique(all_sim$segment_name))
        results = readRDS(paste0("results_simulation",S,".rds"))
        best_K = results$best_K

        Path = single_directory
        MAE_files = list.files(path = Path, pattern = "MAE", recursive = FALSE, full.names = TRUE) 
        RI_files = list.files(path = Path, pattern = "RI", recursive = FALSE, full.names = TRUE) 

        for (j in 1: length(MAE_files)){
            if (j == best_K){
                MAE = readRDS(MAE_files[j])
                RI = readRDS(RI_files[j])

                MAE_list = c(MAE_list, MAE)
                print(MAE_list)
                RI_list = c(RI_list, RI)                
            }  
        }
        setwd("../..")
    }

    avg_MAE = mean(MAE_list)
    print(avg_MAE)
    avg_RI = mean(RI_list)     

    my_list$MAE = c(my_list$MAE, avg_MAE)
    my_list$RI = c(my_list$RI, avg_RI)

    setwd(original_dir)
    print(paste0("directory at the end of the loop sor the first simulation ",getwd()))

}





# my_list <- list(rating=1:4,
#                 animal=c('koala', 'hedgehog', 'sloth', 'panda'),
#                 country=c('Australia', 'Italy', 'Peru', 'China'),
#                 avg_sleep_hours=c(21, 18, 17, 10))
# print(my_list)

# super_sleepers <- data.frame(my_list)
# print(super_sleepers)