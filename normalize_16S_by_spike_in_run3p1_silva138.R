##########################
# GOAL
##########################
# normalize 16S by spike in recovered to estimate "absolute" abundance of each ASV in copies/mL
# normalization methods are based on Lin et al. 2019, Applied and Envrionmental Microbiology
# use taxonomy from silva 138

##########################
# set up environment
##########################
setwd("~/Documents/NYU_drive_local_4_24/P18/Manuscripts/biogeography/Drafts/science/science_revision_1/code")
raw.df <- read.csv(file = "P18_16S_run3p1_ASV_nonplastid_v5_silva138_pre_normalization.csv", header = T)
library(lubridate)
library(ggplot2)
library(reshape2)
library(plyr)
library(stringr)

##########################
# get a column for T. thermophilus reads recovered 
# in each sample (sequencing library)
##########################

# make a dataframe containing just T. thermophilus reads recovered
spike_in.df <- raw.df[grep("s__Thermus_thermophilus", raw.df$Taxon_silva),]
names(spike_in.df)
# sum raw T. thermophilus counts by sequencing library to get total spike in 
# counts recovered in each library
spike_in.tots <- data.frame(tapply(spike_in.df$counts, spike_in.df$libname, sum))
names(spike_in.tots) <- c("spike_in_reads_recovered")
spike_in.tots$libname <- row.names(spike_in.tots) # make a column for the sequencing library name

# merge per-library spike-in totals into raw count spreadsheet
raw.df <- merge(raw.df, spike_in.tots, by="libname")
range(raw.df$spike_in_reads_recovered)

##########################
# specify amount of T. thermophilus DNA spiked into 
# libraries prior to sequencing
##########################
# less DNA was spiked into deeper samples
# due to lower biomass

# initialize a column for mass of T. thermophilus DNA spiked in (ng)
raw.df$ng_spike_in <- NA

# specify spike in amount based on sample depth 
raw.df[which(raw.df$CTDPRS < 200), "ng_spike_in"] <- 0.5 
# samples above 200m received 0.5 ng of spike
raw.df[which(raw.df$CTDPRS > 200 & raw.df$CTDPRS < 500), "ng_spike_in"] <- 0.2 
# samples between 200-500m received 0.2 ng of spike
raw.df[which(raw.df$CTDPRS > 500), "ng_spike_in"] <- 0.1 
# samples below 500m received 0.1 ng of spike

# check that all samples have been assigned spike-in quantities
no_spike <- raw.df[which(is.na(raw.df$ng_spike_in)),]
unique(no_spike$itag_sample_no) 
# only controls are missing spike-in value; good

##########################
# calculate the number of T.therm reads spiked-in for each sample
##########################
# this is based on the ng of DNA spiked-in, 
# the rRNA gene copy number in T. thermophilus, 
# and the length of the T. thermophilus genome
# The full formula is 
# reads spiked in = (ng DNA* 6.022*10^23 * rrna copy number)/ (length of spike-in genome)* (1*10^9 ng/g)* 650g/bp*mol
# use an rRNA copy number of 110; true number varies from 100-120

# proof of concept; check math using values from Lin et al. 2018
# check I get the same number for Cs of T. thermophilus as in original paper:
(14.85*(6.022*10^23) *2)/((2.13*10^6)*(1*10^9)*650) #correct

#now apply to our T. therm spike ins
raw.df$copies_spiked_in <- (raw.df$ng_spike_in*(6.022*10^23) *2)/((2.13*10^6)*(1*10^9)*650)

# check the range of values
range(raw.df$ng_spike_in, na.rm = T)
range(raw.df$copies_spiked_in, na.rm = T)

##########################
# normalize raw counts to T. thermophilus reads
##########################
# estimate abundance of each ASV in gene copies per mL of seawater
# formula is:
# normalized counts ASVi= ((raw reads ASVi in sample)*(reads T. thermophilus spiked in))/((reads T. thermophilus recovered)*(volume of seawater filtered in mL))

# get volume of seawater filtered in mL instead of liters
raw.df$mL_filtered <- raw.df$liters_filtered*1000

# estimate rRNA copies/mL
raw.df$est_copies_per_mL <- (raw.df$counts*raw.df$copies_spiked_in)/(raw.df$spike_in_reads_recovered*raw.df$mL_filtered)

# check range
range(raw.df$est_copies_per_mL, na.rm = T)

# would say "INF" (infinity) wherever there were no reads recovered; change to "0" if applicable
length(which(raw.df$est_copies_per_mL==Inf))
raw.df[which(raw.df$est_copies_per_mL==Inf),"est_copies_per_mL"] <- 0

##########################
# adjust value based on extraction pipeline
##########################

# multiply by 2 because we cut the filter in half
raw.df$est_copies_per_mL <-raw.df$est_copies_per_mL *2

# multiply by 30/4 to account for only 4 uL of the 30 uL total RNA extract being used for cDNA synthesis
raw.df$est_copies_per_mL <- raw.df$est_copies_per_mL *30/4

##########################
# visualize recovery rate across depth
##########################
# make a data frame with just spike-ins recovered
spikes.df <- raw.df[grep("s__Thermus_thermophilus", raw.df$Taxon_silva),]
groupCols <- c("libname", "itag_sample_no", "copies_spiked_in", "CTDPRS")

# summarize by CTD pressure (proxy for depth)
spikes.ply <- ddply(spikes.df, groupCols, function(x) colSums(x[which(names(x)=="counts")], na.rm = T)) 
# get % recovered
spikes.ply$percent_recovered <- spikes.ply$counts/spikes.ply$copies_spiked_in*100

plot1 <- ggplot(spikes.ply, aes(x = percent_recovered, y= -CTDPRS, color = itag_sample_no))+
  geom_point()+
  theme(legend.position = "none")+
  xlab("percent spike-in recovered")+
  theme_bw()+
  ggtitle("Percent spike-in recovered with depth")
plot1
ggsave(filename = "P18_16Sr3_perc_spike_in_recovered.png", plot = plot1, units = "in", width = 4, height = 4)  

# save percent of spike recovered for our records
write.csv(file = "P18_16S_run3p1_perc_spikes_recovered.csv", x= spikes.ply, row.names = F)

##########################
# finalize normalized count sheet
##########################

# remove T.thermophilus reads
raw.df <- raw.df[-grep("s__Thermus_thermophilus", raw.df$Taxon_silva),]

# write out estimated absolute abundance ASV spreadsheet
unique(raw.df[which(is.na(raw.df$est_copies_per_mL)), "itag_sample_no" ]) # only controls are missing values

names(raw.df)

# make output sheet
out.df <- raw.df[,c("libname", "itag_sample_no", "Feature.ID", "Taxon_silva",
                   "Confidence_silva", "est_copies_per_mL")]
# remove controls
out.df <- out.df[complete.cases(out.df),]

# write out
write.csv(file = "P18_16S_run3p1_tax_ASV_nonplastid_silva138_est_cpml.csv", x= raw.df, row.names = F)
