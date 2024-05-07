#1.--------------------------- INSTALLATIONS/DATA IMPORT/DATA PRE-PROCESSING--------------#

#install.packages("Fragman", type = "source")
library(Fragman) ##----loads package


#set up working directory
setwd("C:/Users/kshresth/MSI_analysis/plates")

#My TAILOR-MADE FUNCTIONS for mutant scoring
source("C:/Users/kshresth/MSI_analysis/R_functions/score_easy_kul.R") #CALLING MY FUNCTION which will score only highest main fragmant peak
source("C:/Users/kshresth/MSI_analysis/R_functions/score_easy_kul_with_shutters.R") #CALLING MY FUNCTION which will scores main and the stutter
source("C:/Users/kshresth/MSI_analysis/R_functions/score_easy_kul_A27_125bp_to_160.R") #CALLING MY FUNCTION score_easy_kul which will scores stutter aswell, basically just change in folder name
source("C:/Users/kshresth/MSI_analysis/R_functions/score_easy_kul_A27_125bp_to_150_with_shutter.R") #CALLING MY FUNCTION score_easy_kul which will scores stutter aswell, basically just change in folder name


#import data: 
working_files <- c( "C:/Users/kshresth/MSI_analysis/fragment_plates/Plate_001")

#sroting.inds extracts the channel information containing the fluorescent intensities from the DNA capillary electrophoresis for colors.
my.data <- storing.inds(working_files, channels=5, lets.pullup=FALSE,fourier=FALSE, saturated=FALSE)

#Matching the ladder with peaks found using the ladder.info.attach function
 my.ladder <- c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490,500)

ladder.info.attach(stored=my.data, ladder=my.ladder, channel.ladder=5,
   method="iter2", ladd.init.thresh= 50, 
   draw= TRUE)

#2. ------------------------ANALYSES---------------------------------#

#generate electropherogram plots for each sample 
Plot_name <- basename(getwd())

#----------BLUE PANEL-----------------------#
my.panel_blue <- c(150, 149, 148, 147, 146, 145, 144, 143, 142, 141, 140, 139, 138, 137, 136)

#Get the shutter peaks for microsatellite authentication:

Blue_result_shutter <- score_easy_kul_with_shutters(my.inds=my.data, cols = 1, panel=my.panel_blue,ladder=my.ladder, ploidy=8, left.cond=c(0.3,0), right.cond= 0, warn=TRUE, init.thresh = 150)
dev.off()
shutter_sizing_blue <- get.scores(Blue_result_shutter)
write.table(shutter_sizing_blue, row.names = TRUE, col.names = TRUE,sep=";", file = file.path("C:/Users/kshresth/back up C drive_ 29.8.2018 KSHRESTH/MeioMSI/Results/Peak_calling_score_easy", paste( Plot_name,"_","shutter_sizing_A33.csv")))

#mutant scoring for microsatellite A33

Blue_result_mutants <- score.easy_kul(my.inds=my.data, cols = 1, panel=my.panel_blue,ladder=my.ladder, ploidy=1, left.cond=c(0.3,0.5), right.cond=0, warn=TRUE, init.thresh = 150)
dev.off()
mutant_scoring_blue <- get.scores(Blue_result_mutants)
write.table(mutant_scoring_blue, row.names = TRUE, col.names = TRUE,sep=";", file = file.path("C:/Users/kshresth/MSI_analysis/Results/Peak_calling_score_easy", paste( Plot_name,"_","mutant_scoring_A33.csv")))

#prep for exporting files:

mutant_allele_size <- mutant_scoring_blue - 146
mutant_allele_grouping <- ifelse(abs(mutant_allele_size) <= 0.8, c("normal allele"), c("mutant allele"))
list_to_export_blue <- cbind(shutter_sizing_blue,mutant_scoring_blue , mutant_allele_size)
list_to_export_blue <- round( final_list_to_export_blue, digits = 1)
final_list_to_export_blue <- cbind (list_to_export_blue, mutant_allele_grouping)
colnames(final_list_to_export_blue) = c("peak 1","peak 2", "peak 3", "peak 4", "peak 5", "peak 6", "peak 7", "peak 8","microsatellite allele size_A33 in bp", "mutant allele shift in bps_A33", "allele grouping_A33") # giving columns name

#exporting as .csv:

write.table( final_list_to_export_blue, sep=";", col.names=NA,row.names= TRUE, file = file.path("C:/Users/kshresth/MSI_analysis/Results/individual_loci_mutant_scoring", paste( Plot_name,"_","mutant_scoring_A33.csv")))


#----------GREEN PANEL-----------------------#

my.panel_green <- c(125.5,126.5,127.5,128.5,129.5,130.5,131.5,132.5,133.5,134.5,135.5,136.5,137.5,138.5,139.5,140.5,141.5,142.5,143.5,144.5,146.5,147.5,148.5)

#Get the shutter peaks for microsatellite authentication:

Green_result_shutter <- score_easy_kul_with_shutters (my.inds=my.data, cols = 2, panel=my.panel_green,ladder=my.ladder, ploidy=8, left.cond=c(0.3,0.5), right.cond= 0, warn=TRUE, init.thresh = 150)
dev.off()
shutter_sizing_green <- get.scores(Green_result_shutter)


#mutant scoring for microsatellite A37

Green_result_mutants <- score.easy_kul(my.inds=my.data, cols = 2, panel=my.panel_green,ladder=my.ladder, ploidy=1, left.cond=c(0.3,0.5), right.cond=0, warn=TRUE, init.thresh = 150)
dev.off()
mutant_scoring_green <- get.scores(Green_result_mutants)
write.table(mutant_scoring_green, row.names = TRUE, col.names = TRUE,sep=";", file = file.path("C:/Users/kshresth/MSI analysis/Results/Peak_calling_score_easy", paste( Plot_name,"_","mutant_scoring_A27.csv")))

#prep for exporting files:

mutant_allele_size <- mutant_scoring_green - 137.5
mutant_allele_grouping <- ifelse(abs(mutant_allele_size) <= 0.8, 
     c("normal allele"), c("mutant allele"))

list_to_export_green <- cbind(shutter_sizing_green,mutant_scoring_green , mutant_allele_size)
list_to_export_green <- round(list_to_export_green, digits = 1)
final_list_to_export_green <- cbind (list_to_export_green, mutant_allele_grouping)

colnames(final_list_to_export_green) <- c("peak 1" , "peak 2", "peak 3", "peak 4", "peak 5", "peak 6", "peak 7", "peak 8","microsatellite allele size_A27 in bp", "mutant allele shift in bps_A27", "allele grouping_A27") # giving columns name

#exporting as .csv:
write.table( final_list_to_export_green, sep=";", col.names=NA,row.names= TRUE, file = file.path("C:/Users/kshresth/MSI analysis/Results/individual_loci_mutant_scoring", paste( Plot_name,"_","mutant_scoring_A27.csv")))



#----------YELLOW PANEL-----------------------#


my.panel_yellow <- c(138,140,142,144,145,146,147,148,150)

#Get the shutter peaks for microsatellite authentication:

yellow_result_shutter <- score_easy_kul_with_shutters (my.inds=my.data, cols = 3, panel=my.panel_yellow,ladder=my.ladder, ploidy=2, left.cond=c(0.2, 1.5), right.cond= 0.5, warn=TRUE, init.thresh =500, shift = 1.5)
dev.off()
shutter_sizing_yellow <- get.scores(yellow_result_shutter)



#mutant scoring for microsatellite D14Mit15

yellow_result_mutants <- score.easy_kul (my.inds=my.data, cols = 3, panel=my.panel_yellow,ladder=my.ladder, ploidy=1, left.cond=c(0.2,1.5), right.cond=0.5, warn=TRUE, init.thresh = 500, shift= 1.5)
dev.off()
mutant_scoring_yellow <- get.scores(yellow_result_mutants)

#prep for exporting files:

mutant_allele_size <- mutant_scoring_yellow - 148
mutant_allele_grouping <- ifelse(abs(mutant_allele_size) <= 1.5,c("normal allele"), c("mutant allele"))

list_to_export_yellow = cbind(shutter_sizing_yellow,mutant_scoring_yellow , mutant_allele_size)
list_to_export_yellow <- round(list_to_export_yellow, digits = 1)
final_list_to_export_yellow = cbind (list_to_export_yellow, mutant_allele_grouping)

colnames( final_list_to_export_yellow) = c("peak 1" , "peak 2","microsatellite allele size_D14Mit15 in bp", "mutant allele shift in bps_D14Mit15", "allele grouping_D14Mit15") # giving columns name

write.table( final_list_to_export_yellow, sep=";", col.names=NA,row.names= TRUE, file = file.path("C:/Users/kshresth/back up C drive_ 29.8.2018 KSHRESTH/MeioMSI/Results/individual_loci_mutant_scoring", paste( Plot_name,"_","mutant_scoring_D14Mit15.csv")))

#--------------Quality control rechecking the ladder for every sample------------------#

my.panel_ladder<-c(138,138.5,139,139.5,140,140.5,141,148.5,149,149.5,150,150.5,158.5,159,159.5,160,160.5)
my_ladder_check <- score.easy_kul (my.inds=my.data, cols = 5, panel=my.panel_ladder,ladder=my.ladder, electro=FALSE, init.thresh=100, ploidy=3, left.cond=c(0.3,0), right.cond=0)
dev.off()
score_ladder <- get.scores(my_ladder_check)
colnames(score_ladder) = c("Ladder 139" , "Ladder 150", "Ladder 160") # giving columns name
write.table(score_ladder, row.names = TRUE, col.names = NA,sep=";", file = file.path("C:/Users/kshresth/MeioMSI/Results", paste( Plot_name,"score_ladder.csv")))


#3.---FINAL REPORT, SUMMANRY TABLE OF ALL THE DYES AND LADDER SCORES ----------#


mutation_score_labber_combind = cbind( final_list_to_export_blue, final_list_to_export_green, final_list_to_export_yellow,score_ladder)
write.table(mutation_score_labber_combind, row.names = TRUE, col.names = NA ,sep=";", file = file.path("C:/Users/kshresth/MSI_analysis/Results/Combine_mutant_scoring_for_all_dyes_with_ladder", paste( Plot_name,"_","mutation_score_ladder_combind .csv")))

 



