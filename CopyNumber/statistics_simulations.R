files  = (00,01,02,03)


dir_names = paste0("../",1:12)
names = paste0(1:12)

my_list <- list(simulation = names,
                    MAE = c(),
                    RI = c()
                    )


dir_original=getwd()
# for loop on simulations
for (i in dir_names) {

    setwd(i)
    print(i)
    dir_number=getwd()
    print(dir_number)
    # for loop on iteration of simulation

        MAE_list = c ()
        RI_list = c ()

    for (i in 1:3){ #20

        setwd(paste0("simulation_iteration_",i,"/results"))
        dir = getwd()
        # sapply(BiblioDir,function(dir){length(list.files(dir,pattern='pdf'))})


        all_sim = readRDS("all_sim_input_prepare_input_data.rds")
        S = length(unique(all_sim$segment_name))
        results = readRDS("results_simulation",S,".rds")
        best_K = results$best_K

        BiblioPath = dir
        # BiblioDir = list.dirs(path = BiblioPath, full.names = TRUE, recursive = FALSE)
        MAE_files = list.files(path = BiblioPath, pattern = "MAE", recursive = FALSE, full.names = TRUE) 
        RI_files = list.files(path = BiblioPath, pattern = "RI", recursive = FALSE, full.names = TRUE) 

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

    setwd(dir_original)

}





# my_list <- list(rating=1:4,
#                 animal=c('koala', 'hedgehog', 'sloth', 'panda'),
#                 country=c('Australia', 'Italy', 'Peru', 'China'),
#                 avg_sleep_hours=c(21, 18, 17, 10))
# print(my_list)

# super_sleepers <- data.frame(my_list)
# print(super_sleepers)