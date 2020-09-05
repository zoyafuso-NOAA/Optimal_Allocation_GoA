#########################
## Knitting together optimization results
## If method == 'spatiotemporal'
#########################
rm(list = ls())

############################
## Set up directories
#############################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[2]
VAST_model <- "6g"
github_dir <- paste0(c('/Users/zackoyafuso/Documents/', 
                       'C:/Users/Zack Oyafuso/Documents/',
                       'C:/Users/zack.oyafuso/Work/')[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

load(paste0(github_dir, 'optimization_data.RData'))

####################
## Data objects
####################
master_res_df <- data.frame(id = 1:N)
master_settings <- data.frame()
master_strata_list <- master_strata_stats_list <- list()

for(ifile in dir(paste0(github_dir, 'Spatiotemporal_Optimization/'), 
                 pattern = 'strata', full.names = TRUE)){
        
        #Load object
        load(ifile)
        
        istrata = gsub(ifile, pattern = paste0(github_dir,
                                               'Spatiotemporal_Optimization/',
                                               'optimization_'), 
                       replacement = '')
        istrata = gsub(istrata, pattern = '_strata.RData', replacement = '')
        istrata = as.integer(istrata)
        
        if(ncol(res_df) == 2) res_df = data.frame(res_df[,2])
        if(ncol(res_df) > 2) res_df = res_df[,-1]
        
        names(res_df) = NULL
        
        master_res_df = cbind(master_res_df, res_df)
        master_settings = rbind(master_settings, settings)
        
        temp_strata_list = list()
        for(icol in 1:ncol(res_df)){
                temp_strata_list[[icol]] = data.frame(strata_list[1:9 + 9*(icol -1)])
        }
        
        master_strata_list = c(master_strata_list, temp_strata_list)
}

settings = master_settings[order(master_settings$nstrata),]
res_df = master_res_df[,c(1,1+order(master_settings$nstrata))]
strata_list = master_strata_list[order(master_settings$nstrata)]

strata_to_save = which(settings$nstrata %in% c(5, 10, 15, 20, 30, 60))
settings = settings[strata_to_save,]
names(settings)[2] = 'strata'
settings$id = 1:nrow(settings)

res_df = res_df[, c(1, 1+strata_to_save)]
strata_list = strata_list[strata_to_save]

save(list = c('res_df', 'settings', 'strata_list'),
     file = paste0(github_dir, 'Spatiotemporal_Optimization/',
                   'optimization_knitted_results.RData'))

