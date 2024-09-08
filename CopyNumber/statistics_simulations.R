library(dplyr)
library(formattable)

# files  = (00,01,02,03)
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
dir_names = paste0("../",1:12)
names = paste0(1:12)
# list of average score (among the 20 repetitions of the same simulation type) 
# of the model choosen as best with model selection
my_list <- list(simulation = names,
                    MAE = c(),
                    RI = c(),
                    LOO_K_best_VS_K_real = c(),
                    LogLik_K_best_VS_K_real = c(),
                    Is_K_best_equal_K_real = c()
                    )

original_dir=getwd()
# for loop on simulations
for (i in dir_names) {

    setwd(i)
    print(i)
    simulation_directory=getwd()
    
    # list where to store the values of the best model MAe of a simulation 
    # to then perform the average of these values to put in the main list MAE
    MAE_list = c ()
    RI_list = c ()

    LOO_list = c()
    LogLik_list = c()

    Is_K_best_equal_K_real_list = c()

    simulation_directories = list.dirs(path = simulation_directory, recursive = FALSE, full.names = TRUE) 
    directories_number = length(simulation_directories)

    # for loop on iteration of simulation
    for (i in 1:directories_number){ 

        setwd(paste0("simulation_iteration_",i,"/results"))
        single_directory = getwd()
        # sapply(BiblioDir,function(dir){length(list.files(dir,pattern='pdf'))})

        all_sim = readRDS("all_sim_input_prepare_input_data.rds")
        S = length(unique(all_sim$segment_name))
        results = readRDS(paste0("results_simulation",S,".rds"))
        best_K = results$best_K
        model_selection = readRDS("model_selection.rds")

        LOO_K_best = model_selection %>% filter(K==best_K) %>% pull(LOO)
        LOO_K_real = model_selection %>% filter(K==length(unique(all_sim$tau))) %>% pull(LOO)
        LogLik_K_best = model_selection %>% filter(K==best_K) %>% pull(Log_lik)
        LogLik_K_real = model_selection %>% filter(K==length(unique(all_sim$tau))) %>% pull(Log_lik)        

        Path = single_directory
        MAE_files = list.files(path = Path, pattern = "MAE", recursive = FALSE, full.names = TRUE) 
        RI_files = list.files(path = Path, pattern = "RI", recursive = FALSE, full.names = TRUE) 

        for (j in 1: length(MAE_files)){
            if (j == best_K){
                MAE = readRDS(MAE_files[j])
                RI = readRDS(RI_files[j])


                MAE_list = c(MAE_list, MAE)
                RI_list = c(RI_list, RI)                
            }  
        }

        LOO_list = c(LOO_list, abs(LOO_K_real-LOO_K_best))
        LogLik_list = c(LogLik_list, abs(LogLik_K_real-LogLik_K_best))
        if (best_K == length(unique(all_sim$tau))){
            Is_K_best_equal_K_real_list = c(Is_K_best_equal_K_real_list, TRUE)
        }else{
            Is_K_best_equal_K_real_list = c(Is_K_best_equal_K_real_list, FALSE)       
        }
        setwd("../..")
    }

    avg_LOO_difference = mean(LOO_list)
    avg_LogLik_difference = mean(LogLik_list)
    my_list$LOO_K_best_VS_K_real = c(my_list$LOO_K_best_VS_K_real, avg_LOO_difference)
    my_list$LogLik_K_best_VS_K_real = c(my_list$LogLik_K_best_VS_K_real, avg_LogLik_difference)

    avg_MAE = mean(MAE_list)
    avg_RI = mean(RI_list)     

    my_list$MAE = c(my_list$MAE, avg_MAE)
    my_list$RI = c(my_list$RI, avg_RI)

    my_list$Is_K_best_equal_K_real = c(my_list$Is_K_best_equal_K_real, Mode(Is_K_best_equal_K_real_list) )

    setwd(original_dir)
}


saveRDS(my_list, "score_list.rds")

score_list = readRDS("score_list.rds")
print(score_list)

df_score_list = data.frame(score_list)

table_formattable <- formattable(df_score_list, list(

  MAE = color_tile("white","orange"),
  RI = formatter("span", style = x ~ ifelse(x > 0.70 ,
    style(color = "green", font.weight = "bold"), NA)),
  area(col = c(LOO_K_best_VS_K_real)) ~ normalize_bar("pink", 0.2),
  area(col = c(LogLik_K_best_VS_K_real)) ~ normalize_bar("gray", 0.2),
  Is_K_best_equal_K_real = formatter("span",
    style = x ~ style(color = ifelse(x, "green", "red")),
    x ~ icontext(ifelse(x, "ok", "remove"), ifelse(x, "TRUE", "FALSE")))
))

saveRDS(table_formattable, "table_scores_formattable.rds")

# my_list <- list(rating=1:4,
#                 animal=c('koala', 'hedgehog', 'sloth', 'panda'),
#                 country=c('Australia', 'Italy', 'Peru', 'China'),
#                 avg_sleep_hours=c(21, 18, 17, 10))
# print(my_list)

# super_sleepers <- data.frame(my_list)
# print(super_sleepers)