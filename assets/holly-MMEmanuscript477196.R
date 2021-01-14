
#' The following script generates the figures presented in:
#' 
#' Olsen, E., Eide, C.H., Nilsen, I., Perryman, H.A. and Vikebø, F., 2019. Ecological effects and ecosystem 
#' shifts caused by mass mortality events on early life stages of fish. Frontiers in Marine Science, 6, p.669.

#' 
#' This script is to share the base code for producing the graphs presented in the paper.
#'

#' 
#' The settings and corresponding play data set are just for 
#' executing the code - the resutls do not correspond to the paper. 
#'  





#'
#' Initialization
#' 

#' load libraries   

library(tidyverse)
library(RColorBrewer)
library(extrafont)
#
library(ncdf4)
#loadfonts(device="win") 



# this function was made by Holly based on discussions with Ina and Cecilie
get_PrimaryProduction <- function(fids, 
                                  PRODnc, 
                                  BGMfile){
  nc <- nc_open(PRODnc)
  bboxes <- get_boundary(boxinfo = load_box(bgm = file.path(BGMfile)))
  boxdetails <- bgmfile(BGMfile)
  forout <- matrix(NA, 
                   nrow = dim(ncvar_get(nc,paste(as.character(nordic.group.data$Name[fids[1]]),"Prodn", sep ="")))[2], 
                   ncol = length(fids))
  for (spp in c(1:length(fids))) {
    temp <- ncvar_get(nc,paste(as.character(nordic.group.data$Name[fids[spp]]),"Prodn", sep =""))
    temp=temp*5.7*20/(1000*1000*1000) # Convert from (mg N m-3 d-1) to (tons m-3 d-1) 
    prod.tot <- matrix(NA,nrow=dim(temp)[1],ncol=dim(temp)[2])
    for(i in c(1:dim(temp)[1])){ # Convert from (tons m-3 d-1) to (tons d-1) 
      prod.tot[i,]=temp[i,] * as.numeric(boxdetails$boxes[i,4]) * -as.numeric(boxdetails$boxes[i,3])
    } 
    temp <- temp[!(c(1:dim(temp)[1]) %in% (bboxes + 1)),] 
    forout[,spp] <- colSums(temp)
  }
  PPout <- rowSums(forout)
  return(PPout)
}


ExtractFromNCatltools <- function(fids, # path to fid csv file
                                  init, # path to initialization nc file 
                                  prm_biol, # path to biological prm file
                                  plgn, # path to atlantis input for bgm map
                                  nc_out, # path to output nc file
                                  prm_run, # path to run prm file
                                  fids_age, # selected age-structure fids
                                  fids_pool){ # selected biomass pool fids
  print("!!!!!")
  print("Please double check the bat file to make sure all of the provided paths correspond to those used in the simulation.")
  print("!!!!!")
  # get benthic groups from input nc file
  bps <- load_bps(fgs = file.path(fids), init = file.path(init))
  # get bio_conv from biological prm file
  bio_conv <- get_conv_mgnbiot(prm_biol = file.path(prm_biol))
  # get boundary boxes from bgm file
  bboxes <- get_boundary(boxinfo = load_box(bgm = file.path(plgn)))
  # 1 - Get nc data to calc biomass for age structured groups:
  # Get Nums from nc file
  Nums <- load_nc(nc = file.path(nc_out),
                  bps = bps,
                  fgs = file.path(fids),
                  prm_run = file.path(prm_run),
                  bboxes = bboxes,
                  select_groups = fids_age,
                  select_variable = "Nums")
  # get StructN from nc file
  sn  <- load_nc(nc = file.path(nc_out),
                 bps = bps,
                 fgs = file.path(fids),
                 prm_run = file.path(prm_run),
                 bboxes = bboxes,
                 select_groups = fids_age,
                 select_variable = "StructN")
  # get ResN from nc file
  rn  <- load_nc(nc = file.path(nc_out),
                 bps = bps,
                 fgs = file.path(fids),
                 prm_run = file.path(prm_run),
                 bboxes = bboxes,
                 select_groups = fids_age,
                 select_variable = "ResN")
  # 2 - Get nc data to calc biomass for biomass pooled groups:
  # get N from nc file
  n   <- load_nc(nc = file.path(nc_out), 
                 bps = bps, 
                 fgs = file.path(fids),
                 prm_run = file.path(prm_run), 
                 bboxes = bboxes,
                 select_groups = fids_pool, 
                 select_variable = "N")
  # get physics, specifically vol and dz, from nc file
  vol <- load_nc_physics(nc = file.path(nc_out),
                         prm_run = file.path(prm_run), 
                         bboxes = bboxes, 
                         aggregate_layers = FALSE,
                         select_physics = c("volume", "dz"))
  # CALC THE BIOMASS OF FIDS ACROSS SPACE
  print("Calculating group(s) biomass using calculate_biomass_spatial from atlantistools")
  B <- calculate_biomass_spatial(nums = Nums, 
                                 sn = sn, 
                                 rn = rn, 
                                 n = n, 
                                 vol_dz = vol, 
                                 bio_conv = bio_conv, 
                                 bps = bps)
  print("Calculating complete.")
  print("Merging output into one dataframe. This may take a couple minutes...")
  list(B,rn,sn,Nums)
}


#' load model inputs  
#' - in the working directory there needs to be a folder with the necessary model inputs

#' list of functional group codes linked to full name 
nordic.group.data <- read.table("input_files/nordic_groups_Oil.csv", dec=".", sep=",", header=TRUE)
nordic.group.data <- nordic.group.data[which(nordic.group.data$IsTurnedOn == 1),] # REMOVE FIDS THAT ARE TURNED OFF

#' polygon data
File_bgm <- "input_files/Nordic02.bgm" 
#bboxes <- get_boundary(boxinfo = load_box(bgm = file.path(File_bgm))) # I AM NOT SURE IF I NEED THIS!!!! 

#' run data
toutinc <- 73 # DECLARING LENGTH BETWEEN TIMESTEPS (SEE THE RUN PRM FILE)
eventyr <- 61 # The time step just before the infliction of the mortality event 

#' basic info of functional groups for computing indicators 
basic=read.table("input_files/NOBA_BasicInfo_Oil.csv",header=T,sep=';',na.strings = "NA",stringsAsFactors = F)

#' load model ouputs  
#' - in the working directory there needs to be a folder with the output files
#' - this directory consists of individual folders for each scenario run
#' - inside each of these folders are the necessary output files. 

#' identify all of the individual output folders within the ouput_files directory 
#' - a subset of the folders is provided to test run the code
folder0  = "output_00_26.01.19"
#folder1  = "output_01_26.01.19"
#folder2  = "output_02_26.01.19"
#folder3  = "output_03_26.01.19"
#folder4  = "output_04_26.01.19"
#folder5  = "output_05_26.01.19"
#folder6  = "output_06_26.01.19"
#folder7  = "output_07_26.01.19"
#folder8  = "output_08_26.01.19"
#folder9  = "output_09_26.01.19"
folder10 = "output_10_26.01.19"
folder11 = "output_11_26.01.19"
folder12 = "output_12_26.01.19"

#' read in all the biomass data
run_00 <- read.table(paste("output_files/", folder0, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
#run_01 <- read.table(paste("output_files/", folder1, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
#run_02 <- read.table(paste("output_files/", folder2, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
#run_03 <- read.table(paste("output_files/", folder3, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
#run_04 <- read.table(paste("output_files/", folder4, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
#run_05 <- read.table(paste("output_files/", folder5, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
#run_06 <- read.table(paste("output_files/", folder6, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
#run_07 <- read.table(paste("output_files/", folder7, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
#run_08 <- read.table(paste("output_files/", folder8, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
#run_09 <- read.table(paste("output_files/", folder9, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_10 <- read.table(paste("output_files/", folder10, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_11 <- read.table(paste("output_files/", folder11, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_12 <- read.table(paste("output_files/", folder12, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)

#' read in all the data from the nc file
nc_00 <- nc_open(paste("output_files/", folder0, "/ina_results_001.nc", sep=""))
#nc_01 <- nc_open(paste("output_files/", folder1, "/ina_results_001.nc", sep=""))
#nc_02 <- nc_open(paste("output_files/", folder2, "/ina_results_001.nc", sep=""))
#nc_03 <- nc_open(paste("output_files/", folder3, "/ina_results_001.nc", sep=""))
#nc_04 <- nc_open(paste("output_files/", folder4, "/ina_results_001.nc", sep=""))
#nc_05 <- nc_open(paste("output_files/", folder5, "/ina_results_001.nc", sep=""))
#nc_06 <- nc_open(paste("output_files/", folder6, "/ina_results_001.nc", sep=""))
#nc_07 <- nc_open(paste("output_files/", folder7, "/ina_results_001.nc", sep=""))
#nc_08 <- nc_open(paste("output_files/", folder8, "/ina_results_001.nc", sep=""))
#nc_09 <- nc_open(paste("output_files/", folder9, "/ina_results_001.nc", sep=""))
nc_10 <- nc_open(paste("output_files/", folder10, "/ina_results_001.nc", sep=""))
nc_11 <- nc_open(paste("output_files/", folder11, "/ina_results_001.nc", sep=""))
nc_12 <- nc_open(paste("output_files/", folder12, "/ina_results_001.nc", sep=""))

# Colors for plots
#http://research.stowers.org/mcm/efg/R/Color/Chart/ColorChart.pdf
cCO10<-"#BEC1D4" ;cCO50<-"#7D87B9" ;cCO90<-"#023FA5" 
cHA10<-"#F0F2D5" ;cHA50<-"#B7D0A0" ;cHA90<-"#57945C" 
cHE10<-"#F1DE81" ;cHE50<-"#EEAB65" ;cHE90<-"#B9534C" 
cMX10<-"#D8BFD8" ;cMX50<-"#BA55D3" ;cMX90<-"#800080" 


#'
#'
#' Figure 2
#' - this is also the code used to make Figure S2
#' 
#'

#' make dataframe for plot
biom.df <-bind_rows(`BASE` = run_00,
                    #`CO10` = run_01,
                    #`CO50` = run_02,
                    #`CO90` = run_03,
                    #`HE10` = run_04,
                    #`HE50` = run_05,
                    #`HE90` = run_06,
                    #`HA10` = run_07,
                    #`HA50` = run_08,
                    #`HA90` = run_09,
                    `MX10` = run_10,
                    `MX50` = run_11,
                    `MX90` = run_12,
                    .id="Runs")%>%
  gather(key=Species_code, value=Biomass, -Time, -Runs) %>%
  rename(Day =Time) %>% 
  mutate(Year=Day/365) %>% 
  mutate(Season = Year %% 1) %>%
  filter(Season==0)  %>%
  filter(Year>eventyr) %>% 
  filter(Species_code %in% c("HAD", "SSH", "NCO")) %>%   # Choose species
  left_join(nordic.group.data, by=c("Species_code"="Code"))         # Add complete name of species

# adjust year so you are ploting years since mass mortality event
biom.df$Year <- biom.df$Year - eventyr

#' set order of runs, species and colours for plot
run.order <- c(#"CO10","CO50","CO90","HA10","HA50","HA90","HE10","HE50","HE90",
              "MX10","MX50","MX90","BASE")
sp.order <- c("Northeast Arctic Cod","Northeast Arctic Haddock","Norwegian Spring Spawning Herring")
colours <- c(#"#a6cee3","#80b1d3","#1f78b4","#b3de69","#33a02c","#006837","#fdb462","#ff7f00","#b15928",
             "#cab2d6","#bc80bd","#6a3d9a","#000000")

#' make plot
plot <- biom.df %>% 
  filter(Species_code %in% c("SSH","HAD","NCO")) %>% 
  mutate(Runs = factor(Runs, levels = run.order)) %>%
  arrange(Runs) %>% 
  mutate(LongName = factor(LongName, levels = sp.order)) %>%
  arrange(LongName)  %>% 
  mutate(Biomass=Biomass/250000) %>%                                   # Create index of biomass
  ggplot(aes(x=Year, y=Biomass, color=Runs)) +  
  geom_line(lwd=1.5) +  
  scale_colour_manual(values=colours)+  
  labs(x="Time (Year post MME)", y="Biomass index", title="", colour="Scenarios") +
  theme_minimal() +
  theme(#text=element_text(family="Calibri",size=25),
        plot.title = element_text(size=24)) +
  facet_wrap(~LongName, scales = "free_y") # Plot species seperatly (different scale on y-axis)
plot

#' save plot
# ggsave("Figure_2.tiff", units="in", height=5.5, width=17, dpi=300 )
dev.off()

#' remove misc
remove(biom.df,plot,run.order,sp.order,colours)




#'
#'
#' Figure 3
#' - there are warnings but the functions still work
#'

#' create timeframe
time=ncvar_get(nc_00,"t")

#' Choose species
sel.Fish=c("North_atl_cod","Haddock","Norwegian_ssh") 

#' make dataframes for runs
#' - this is done for each selected run
#' - the base run is addes as a column in the dataframes to compute values relative to the baseline

#' baseline (run 00)

Biomass_00 <- lapply(sel.Fish, function(f) {          # All age-structured groups have
  if (f == "Capelin") {                               # 10 age classes, except:
    ycl = c(1:5)                                      # Capelin, snow crab and sperm whale
  } else if (f == "Snow_crab") {
    ycl = c(1:6)
  } else if (f == "Sperm_whale") {
    ycl = c(1:8)
  } else{
    ycl = c(1:10)
  }
  sapply(ycl, function(y){ 
    nums = ncvar_get(nc_00, paste(f, y, "_Nums", sep = ""))      # Extract numbers from nc-file
    struct = ncvar_get(nc_00,paste(f, y, "_StructN", sep = ""))  # Extract structural weight
    res = ncvar_get(nc_00,paste(f, y, "_ResN", sep = ""))        # Extract reserved weight
    biom=(struct+res)*nums*5.7*20/1e9                           # Convert from Nmg to tonnes
    biom_ts=colSums(colSums(biom))                              #Sums layers and polygons
    biom_ts
    
  }) %>% as.tibble() %>%                                        
    mutate(Species = f) %>% 
    rowid_to_column("time") %>% 
    mutate(time1=time-1) %>% 
    mutate(day = time1*toutinc, Year = day/365) %>% 
    select(-time, -time1)
  
}) %>% bind_rows() %>% 
  gather(key = Age_class, value = Biomass, -Species, -Year) %>% 
  mutate(Age_class = gsub("V","", Age_class), 
         Age_class = as.integer(Age_class)) %>% 
  filter(Year>=eventyr) %>%      # Filter out years before mortality event
  rename(Biomass00=Biomass) # Rename biomass columne to add to other data frames

#' Run nc_10 

Biomass_10 <- lapply(sel.Fish, function(f) {
  if (f == "Capelin") {
    ycl = c(1:5)
  } else if (f == "Snow_crab") {
    ycl = c(1:6)
  } else if (f == "Sperm_whale") {
    ycl = c(1:8)
  } else{
    ycl = c(1:10)
  }
  sapply(ycl, function(y){
    nums = ncvar_get(nc_10, paste(f, y, "_Nums", sep = ""))
    struct = ncvar_get(nc_10,paste(f, y, "_StructN", sep = ""))
    res = ncvar_get(nc_10,paste(f, y, "_ResN", sep = ""))
    biom=(struct+res)*nums*5.7*20/1e9 #Convert from Nmg to tonnes
    biom_ts=colSums(colSums(biom))    #Sums layers and polygons
    biom_ts
    
  }) %>% as.tibble() %>% 
    mutate(Species = f) %>% 
    rowid_to_column("time") %>% 
    mutate(time1=time-1) %>% 
    mutate(day = time1*toutinc, Year = day/365) %>% 
    select(-time, -time1)
  
}) %>% bind_rows() %>% 
  gather(key = Age_class, value = Biomass, -Species, -Year) %>% 
  mutate(Age_class = gsub("V","", Age_class), 
         Age_class = as.integer(Age_class)) %>% 
  filter(Year>=eventyr) %>%  # Filter out years before mortality event
  left_join(Biomass_00) # Add column from base run

#' Run nc_11 

Biomass_50 <- lapply(sel.Fish, function(f) {
  if (f == "Capelin") {
    ycl = c(1:5)
  } else if (f == "Snow_crab") {
    ycl = c(1:6)
  } else if (f == "Sperm_whale") {
    ycl = c(1:8)
  } else{
    ycl = c(1:10)
  }
  sapply(ycl, function(y){
    nums = ncvar_get(nc_11, paste(f, y, "_Nums", sep = ""))
    struct = ncvar_get(nc_11,paste(f, y, "_StructN", sep = ""))
    res = ncvar_get(nc_11,paste(f, y, "_ResN", sep = ""))
    biom=(struct+res)*nums*5.7*20/1e9 #Convert from Nmg to tonnes
    biom_ts=colSums(colSums(biom))    #Sums layers and polygons
    biom_ts
    
  }) %>% as.tibble() %>% 
    mutate(Species = f) %>% 
    rowid_to_column("time") %>% 
    mutate(time1=time-1) %>% 
    mutate(day = time1*toutinc, Year = day/365) %>% 
    select(-time, -time1)
  
}) %>% bind_rows() %>% 
  gather(key = Age_class, value = Biomass, -Species, -Year) %>% 
  mutate(Age_class = gsub("V","", Age_class), 
         Age_class = as.integer(Age_class)) %>% 
  filter(Year>=eventyr) %>%  # Filter out years before mortality event
  left_join(Biomass_00) # Add column from base run

#' Run nc_11

Biomass_90 <- lapply(sel.Fish, function(f) {
  if (f == "Capelin") {
    ycl = c(1:5)
  } else if (f == "Snow_crab") {
    ycl = c(1:6)
  } else if (f == "Sperm_whale") {
    ycl = c(1:8)
  } else{
    ycl = c(1:10)
  }
  sapply(ycl, function(y){
    nums = ncvar_get(nc_11, paste(f, y, "_Nums", sep = ""))
    struct = ncvar_get(nc_11,paste(f, y, "_StructN", sep = ""))
    res = ncvar_get(nc_11,paste(f, y, "_ResN", sep = ""))
    biom=(struct+res)*nums*5.7*20/1e9 #Convert from Nmg to tonnes
    biom_ts=colSums(colSums(biom))    #Sums layers and polygons
    biom_ts
    
  }) %>% as.tibble() %>% 
    mutate(Species = f) %>% 
    rowid_to_column("time") %>% 
    mutate(time1=time-1) %>% 
    mutate(day = time1*toutinc, Year = day/365) %>% 
    select(-time, -time1)
  
}) %>% bind_rows() %>% 
  gather(key = Age_class, value = Biomass, -Species, -Year) %>% 
  mutate(Age_class = gsub("V","", Age_class), 
         Age_class = as.integer(Age_class)) %>% 
  filter(Year>=eventyr) %>%  # Filter out years before mortality event
  left_join(Biomass_00) # Add column from base run

#' compare runs with baseline

#' make combined badaframe of all runs

Biom_all <- bind_rows(`MX10` = Biomass_10,   # (run with 10 % mortality added)
                      `MX50` = Biomass_50,   # (run with 50 % mortality added)
                      `MX90` = Biomass_90,   # (run with 90 % mortality added)
                      .id="Run")%>%
  mutate(Season = Year %% 1) %>%
  filter(Season==0) %>%      # As data is printed 5 times per year we choose the first timestep
  mutate(Diff=(Biomass-Biomass00)/Biomass00*100) #Difference (%) between base and modified runs

#' add black line representing total change of all age classes

Biom_tot <- Biom_all %>% 
  group_by(Run, Species, Year) %>% 
  summarise(Biomass00_sum=sum(Biomass00),
            Biomass_sum=sum(Biomass)) 

Bar_plot <- left_join(Biom_tot,Biom_all) %>% 
  mutate(Change=Diff*(Biomass00/Biomass00_sum)) 

Line_plot <-Biom_tot %>% 
  mutate(Tot_change=(Biomass_sum-Biomass00_sum)/Biomass00_sum*100)

Bar_plot$Year = Bar_plot$Year - eventyr
Line_plot$Year = Line_plot$Year - eventyr

# make plot

Plot <- Bar_plot %>% filter(Species %in% c( "North_atl_cod")) %>%    # Choose species
  filter( Age_class==1 | Age_class==2 | Age_class==3 | Age_class==4 | Age_class==5 | Age_class==6 |Age_class==7| Age_class==8 | Age_class==9 | Age_class==10) %>% 
  ggplot(aes(x=Year, y=Change, fill=as.factor(Age_class))) +
  scale_fill_brewer(palette = "Paired")+  
  geom_col() +
  geom_hline(yintercept = 0, size=1) + # Add base line at 0
  geom_line(aes(x=Year, y=Tot_change),Line_plot %>%   # Add line to plot
              filter(Species %in% c( "North_atl_cod")), inherit.aes = FALSE , size =1) +
  labs(x="Time (Year post MME)", y="Change in biomass relative \n to baseline (%)", title="Northeast Arctic Cod", fill="Age class") +
  theme_minimal()+
  theme(#text=element_text(family="Calibri",size=25), 
        plot.title = element_text(size=30,face="bold")) +
  facet_wrap(~Run)    # Plot all runs seperatly  (same scale on axis)      
Plot 

# save
#ggsave("Figure_3.tiff", units="in", height=6, width=16, dpi=300 )
dev.off()

# remove misc.
remove(Bar_plot,Biom_all,Biom_tot,Biomass_00,Biomass_10,Biomass_50,Biomass_90,Line_plot,time,sel.Fish,Plot)




#' The following also needs these libraries
#' they are not loaded in the begining because
#' there is a conflict

library(ReactiveAtlantis)
library(atlantistools) # need for get_boundary()
library(rbgm) # need for bgmfile()
library(fmsb) # for radar plot




#'
#'
#' Figure 4
#' 
#'

#' Figures 4 and S4 were created by R-script developed by Kelli Johnson, Isaac Kaplan, and Gavin Fay for Olsen et al. (2018). 
#' This code has been published: https://github.com/r4atlantis/common_scenarios_analysis/tree/master/sept17

#' Olsen, E., Kaplan, I. C., Ainsworth, C., Fay, G., Gaichas, S., Gamble, R., Girardin, R., Eide, C. H., Ihde, T. F.,
#' Morzaria-Luna, H. N., Johnson, K. F., Savina-Rolland, M., Townsend, H., Weijerman, M., Fulton, E. A. and Link, 
#' J. S. (2018). Ocean Futures Under Ocean Acidification, Marine Protection, and Changing Fishing Pressures Explored 
#' Using a Worldwide Suite of Ecosystem Models. Frontiers in Marine Science, 5: 64.




#'
#'
#' Figure 5
#' 
#'

#' create list of the indicators for calculations 
indicators_e <- c("PelB.PP", # ECOLOGICAL INDICATORS
                  "Bio.PP",
                  "MTL.B",
                  "PropPred",
                  "DemT.PelT",
                  "Dem.Pel",
                  "Dem.PP",
                  #
                  "TotPC", # FISHERIES INDICATORS 
                  "TotC",
                  "MTL.C",
                  "FishExRt",
                  "ExRt",
                  "Val",
                  "TotFC",
                  "TotDC")

#'
#' RADAR PLOTS 
#' AVERAGES I) 5 YEARS, II) 10 YEARS, AND III) 20 YEARS POST MME
#'
t.span.out <- 5 # Define the temporal span of data going into calculations; 5, 10, or 20

#' seasonality flag (T/F):
#' (T) - you are collecting all of the seasonal data reported in the out files 
#' (F) - you are collecting only the data reported at the end of the year (i.e., every 365 time steps) 
t.season <- F # !!! NOT PROGRAMMED FOR t.season <- T YET !!! 

#' run names
# tier.t=c("CO10","CO50","CO90",
#          "HE10","HE50","HE90",
#          "HA10","HA50","HA90",
#          "MX10","MX50","MX90","BASE")
tier.t=c("MX10","MX50","MX90","BASE")

#' create a list of the out folders to loop through
tier <- c(#folder1, folder2, folder3, folder4, folder5, folder6, folder7, folder8, folder9,
          folder10, folder11, folder12,
          folder0)

#' create array for collecting indicator computations for radar plot
forplot_radarc=array(0,dim=c((length(tier)+2),length(indicators_e)))

# Loop through each model scenario out to collect data for radar plot
# (this method works but is inefficient)
for(i in 1:length(tier)){ 
  biom=read.table(paste("output_files/",tier[i],'/ina_results_001BiomIndx.txt',sep=''),header=T) 
  catch=read.table(paste("output_files/",tier[i],'/ina_results_001Catch.txt',sep=''),header=T) 
  last.yr <- floor(biom$Time[length(biom$Time)]/365) 
  if(t.season){
    # !!! NOT PROGRAMED !!! 
  } else {
    # !!!! the above line is assuming you are collecting data from the last (t.span.out years) of the simulation 
    data.grab <- c((eventyr + 1):(eventyr + t.span.out))
    # !!!! the above line is assuming you are collecting data immediately following the mass mortality event 
  }
  b.data.for.calc <- biom[which(biom$Time %in% (data.grab*365)),] # you need to correct data.grab so it is off the same scale as the txt out files
  c.data.for.calc <- catch[which(catch$Time %in% (data.grab*365)),] # you need to correct data.grab so it is off the same scale as the txt out files
  for(j in 1:length(indicators_e)){
    if(indicators_e[j] == "PelB.PP"){
      guild=which(basic$Atlantis.species.code=='IsPelagic')  # row number in basic that corresponds to IsPelagic 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPelagic
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(b.data.for.calc[,sp])                 # total biomass
      guild=which(basic$Atlantis.species.code=='IsPrimaryProducer')  # row number in basic that corresponds to IsPrimaryProducer 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPrimaryProducer
      sp=sp-1                                                # correct functional group identifer 
      PPout <- get_PrimaryProduction(sp, # NOTE !!!! this function was developed under the assumption that nordic.group.data is in the Environment
                                     paste("output_files/",tier[i],"/ina_results_001PROD.nc",sep=""), 
                                     File_bgm)
      tempPP <- PPout[(biom$Time / 365) %in% (data.grab)]
      tempvalue <- biom_var1 / tempPP
      remove(biom_var1,tempPP,PPout,sp,guild)
    } else if(indicators_e[j] == "Bio.PP"){ 
      biom_var1=rowSums(b.data.for.calc[,c(2:dim(basic)[2])]) # total biomass
      guild=which(basic$Atlantis.species.code=='IsPrimaryProducer')  # row number in basic that corresponds to IsPrimaryProducer 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPrimaryProducer
      sp=sp-1                                                # correct functional group identifer 
      PPout <- get_PrimaryProduction(sp, # NOTE !!!! this function was developed under the assumption that nordic.group.data is in the Environment
                                     paste("output_files/",tier[i],"/ina_results_001PROD.nc",sep=""), 
                                     File_bgm)
      tempPP <- PPout[(biom$Time / 365) %in% (data.grab)]
      tempvalue <- biom_var1 / tempPP
      remove(biom_var1,tempPP,PPout,sp,guild)
    } else if(indicators_e[j] == "MTL.B"){ 
      b.data.for.calc$MTL <- NA
      for (tempyr in c(1:nrow(b.data.for.calc))) {
        TLs <- as.numeric(as.character(basic[which(basic$Atlantis.species.code=='TrophicLevel'),-1]))
        biom.var <- b.data.for.calc[tempyr,c(2:(ncol(basic)-1))]
        numtr <- biom.var * TLs
        b.data.for.calc[tempyr,"MTL"] <- sum(numtr) / sum(biom.var)
      }
      tempvalue <- b.data.for.calc$MTL
      remove(TLs, biom.var, numtr, tempyr)
    } else if(indicators_e[j] == "PropPred"){ 
      guild=which(basic$Atlantis.species.code=='IsPredatoryFish') # row number in basic that corresponds to IsPredatoryFish 
      sp=which(basic[guild,]==1)                                  # functional groups (columns) that meet IsPredatoryFish
      sp=sp-1                                                     # correct functional group identifer 
      biom_var1=rowSums(b.data.for.calc[,sp])                     # total biomass
      guild=which(basic$Atlantis.species.code=='IsFish')     # row number in basic that corresponds to IsFish 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsFish
      sp=sp-1                                                # correct functional group identifer 
      biom_var2=rowSums(b.data.for.calc[,sp])                 # total biomass
      tempvalue <- biom_var1 / biom_var2
      remove(guild,sp,biom_var1,biom_var2)
    } else if(indicators_e[j] == "DemT.PelT"){ 
      guild=which(basic$Atlantis.species.code=='IsDemersalTeleost') # row number in basic that corresponds to IsDemersalTeleost 
      sp=which(basic[guild,]==1)                                    # functional groups (columns) that meet IsDemersalTeleost
      sp=sp-1                                                       # correct functional group identifer 
      biom_var1=rowSums(b.data.for.calc[,sp])                       # total biomass
      guild=which(basic$Atlantis.species.code=='IsPelagicTeleost') # row number in basic that corresponds to IsPelagicTeleost 
      sp=which(basic[guild,]==1)                                   # functional groups (columns) that meet IsPelagicTeleost
      sp=sp-1                                                      # correct functional group identifer 
      biom_var2=rowSums(b.data.for.calc[,sp])                      # total biomass
      tempvalue <- biom_var1 / biom_var2
      remove(guild,sp,biom_var1,biom_var2)
    } else if(indicators_e[j] == "Dem.Pel"){ 
      guild=which(basic$Atlantis.species.code=='IsDemersal') # row number in basic that corresponds to IsDemersalTeleost 
      sp=which(basic[guild,]==1)                                    # functional groups (columns) that meet IsDemersalTeleost
      sp=sp-1                                                       # correct functional group identifer 
      biom_var1=rowSums(b.data.for.calc[,sp])                       # total biomass
      guild=which(basic$Atlantis.species.code=='IsPelagic') # row number in basic that corresponds to IsPelagicTeleost 
      sp=which(basic[guild,]==1)                                   # functional groups (columns) that meet IsPelagicTeleost
      sp=sp-1                                                      # correct functional group identifer 
      biom_var2=rowSums(b.data.for.calc[,sp])                      # total biomass
      tempvalue <- biom_var1 / biom_var2
      remove(guild,sp,biom_var1,biom_var2)
    } else if(indicators_e[j] == "Dem.PP"){ 
      guild=which(basic$Atlantis.species.code=='IsDemersal')  # row number in basic that corresponds to IsDemersal 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsDemersal
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(b.data.for.calc[,sp])                 # total biomass
      guild=which(basic$Atlantis.species.code=='IsPrimaryProducer')  # row number in basic that corresponds to IsPrimaryProducer 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPrimaryProducer
      sp=sp-1                                                # correct functional group identifer 
      PPout <- get_PrimaryProduction(sp, # NOTE !!!! this function was developed under the assumption that nordic.group.data is in the Environment
                                     paste("output_files/",tier[i],"/ina_results_001PROD.nc",sep=""), 
                                     File_bgm)
      tempPP <- PPout[(biom$Time / 365) %in% (data.grab)]
      tempvalue <- biom_var1 / tempPP
      remove(biom_var1,tempPP,PPout,sp,guild)
    } else if(indicators_e[j] == "TotPC"){ 
      guild=which(basic$Atlantis.species.code=='IsPelagic')  # row number in basic that corresponds to IsPelagic 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPelagic
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(c.data.for.calc[,sp])                 # total biomass
      tempvalue <- biom_var1
      remove(biom_var1,sp,guild)
    } else if(indicators_e[j] == "TotC"){ 
      biom_var1=rowSums(c.data.for.calc[,c(2:dim(basic)[2])]) # total biomass
      tempvalue <- biom_var1
      remove(biom_var1)
    } else if(indicators_e[j] == "MTL.C"){ 
      c.data.for.calc$MTL <- NA
      for (tempyr in c(1:nrow(c.data.for.calc))) {
        TLs <- as.numeric(as.character(basic[which(basic$Atlantis.species.code=='TrophicLevel'),-1]))
        biom.var <- c.data.for.calc[tempyr,c(2:(ncol(basic)-1))]
        numtr <- biom.var * TLs
        c.data.for.calc[tempyr,"MTL"] <- sum(numtr) / sum(biom.var)
      }
      tempvalue <- c.data.for.calc$MTL
      remove(TLs, biom.var, numtr)
    } else if(indicators_e[j] == "FishExRt"){ 
      guild1=which(basic$Atlantis.species.code=='IsTarget') # row number in basic that corresponds to IsTarget 
      guild2=which(basic$Atlantis.species.code=='IsFish')   # row number in basic that corresponds to IsTarget 
      sp1=which(basic[guild1,]==1)                          # functional groups (columns) 
      sp2=which(basic[guild2,]==1)                          # functional groups (columns) 
      sp <- sp1[sp1 %in% sp2]                               # targeted groups that are also fish
      sp=sp-1                                               # correct functional group identifer 
      biom_var1=rowSums(c.data.for.calc[,sp])               # total catch
      biom_var2=rowSums(b.data.for.calc[,sp])               # total biomass
      tempvalue <- biom_var1 / biom_var2
      remove(guild1, guild2, sp1, sp2, biom_var1,  biom_var2)
    } else if(indicators_e[j] == "ExRt"){ 
      guild=which(basic$Atlantis.species.code=='IsTarget') # row number in basic that corresponds to IsTarget 
      sp=which(basic[guild,]==1)                             # functional groups (columns) 
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(c.data.for.calc[,sp])                 # total catch
      biom_var2=rowSums(b.data.for.calc[,sp])                 # total biomass
      tempvalue <- biom_var1 / biom_var2
      remove(guild, sp, biom_var1, biom_var2)
    } else if(indicators_e[j] == "Val"){ 
      guild=which(basic$Atlantis.species.code=='USDollarsPerTon') # row number in basic that corresponds to USDollarsPerTon 
      sp=which(is.na(basic[guild,]) == F)                         # functional groups (columns) that that have USDollarsPerTon
      biom_prc <- basic[guild,sp[-1]]                            # prices
      sp=sp[-1]-1                                                 # correct functional group identifer 
      biom_totc=c.data.for.calc[,sp]                              # catch (ton) of fids 
      valprod <- t(t(as.matrix(biom_totc)) * as.numeric(biom_prc))
      biom_var1=rowSums(valprod)                    
      tempvalue <- biom_var1
      remove(guild, sp, biom_prc, biom_totc, valprod, biom_var1)
    } else if(indicators_e[j] == "TotFC"){ 
      guild=which(basic$Atlantis.species.code=='IsFish')  # row number in basic that corresponds to IsFish 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsFish
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(c.data.for.calc[,sp])                 # total biomass
      tempvalue <- biom_var1
      remove(biom_var1,sp,guild)
    } else if(indicators_e[j] == "TotDC"){ 
      guild=which(basic$Atlantis.species.code=='IsDemersal')  # row number in basic that corresponds to IsDemersal 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsDemersal
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(c.data.for.calc[,sp])                 # total biomass
      tempvalue <- biom_var1
      remove(biom_var1,sp,guild)
    } 
    forplot_radarc[(i+2),j]=mean(tempvalue)
  } # end of j loop
  remove(biom,catch,last.yr,b.data.for.calc,c.data.for.calc,tempvalue,data.grab,tempyr)
}

#' Prep data for plot - compute relative values 
forplot <- t(t(forplot_radarc)/forplot_radarc[dim(forplot_radarc)[1],]) 

#' create names - !!! THIS NEEDS TO BE BASED ON THE ORDER IN forplot
colnames(forplot) <- c('Pel.bio/PP',
                       'Bio/PP',
                       'MTL.bio',
                       'Predfish.Prop',
                       'Dem/pel.fish',
                       'Dem/Pel',
                       'Dem.bio/PP',
                       #
                       "Pel.C",
                       "Tot.C",
                       "MTL.C",
                       "Fish.exp.rt",
                       "Exp.rt",
                       "Val",
                       "Fish.C",
                       "Dem.C")
rownames(forplot)=c('max','min',tier.t)

# Create data frame for the radar plots 
forplot1 <- forplot[,c(1:7)]
forplot2 <- forplot[,c(8:15)]

# add min/max values for the first two rows
forplot1[1,] <- 1 #max(forplot,na.rm = T)
forplot1[2,] <- 0.2 #min(forplot[c(1,3:(length(tier)+2)),] ,na.rm = T)
forplot2[1,] <- 1.1 #max(forplot,na.rm = T)
forplot2[2,] <- 0.7 #min(forplot[c(1,3:(length(tier)+2)),] ,na.rm = T)
seqforplot1 <- seq(0.2,1,0.2)
seqforplot2 <- seq(0.7,1.1,0.1)
forplot1=data.frame(forplot1)
forplot2=data.frame(forplot2)

# MAKE PLOTS
#tiff(paste("Radar",t.span.out,".tiff",sep = ""), width = 12, height = 4, units = "in", res = 300)
par(mfrow=c(1,2),xpd=TRUE)
par(#oma=c(0,0,0,0),
  mar=c(0,0,0,6),xpd=TRUE)
# Create radar plot of ecological indicators 
radarchart(forplot1,
           plwd=rep(2,13),
           plty=c(rep(2,4)),
           #plty=c(rep(1,9),rep(2,4)),
           pcol=c(cMX10,cMX50,cMX90,"grey","black"),
           #pcol=c(cCO10,cCO50,cCO90,cHE10,cHE50,cHE90,cHA10,cHA50,cHA90,cMX10,cMX50,cMX90,"grey","black"),
           axistype=1,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seqforplot1, cglwd=0.8,
           vlcex=0.8)
par(mar=c(0,0,0,6),xpd=TRUE)
# Create radar plot of fisheries indicators 
radarchart(forplot2,
           plwd=rep(2,13),
           plty=c(rep(2,4)),
           #plty=c(rep(1,9),rep(2,4)),
           pcol=c(cMX10,cMX50,cMX90,"grey","black"),
           #pcol=c(cCO10,cCO50,cCO90,cHE10,cHE50,cHE90,cHA10,cHA50,cHA90,cMX10,cMX50,cMX90,"grey","black"),
           axistype=1,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seqforplot2, cglwd=0.8,
           vlcex=0.8)
legend("topright",
       #lty=c(rep(1,9),rep(2,4)),
       lty=c(rep(2,4)),
       lwd = rep(2,13),
       #col=c(cCO10,cCO50,cCO90,cHE10,cHE50,cHE90,cHA10,cHA50,cHA90,cMX10,cMX50,cMX90,"grey","black"),
       col=c(cMX10,cMX50,cMX90,"grey","black"),
       legend=tier.t,
       bty = "n",
       inset=c(-0.22,0))
dev.off()

# remove misc
remove(indicators_e,i,j,t.season,t.span.out,tier.t,tier,forplot_radarc,seqforplot1,seqforplot2,forplot,forplot1,forplot2)
