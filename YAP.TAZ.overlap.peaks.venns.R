#needed to install GenomicRanges

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")  # Install BiocManager if not already installed
}
BiocManager::install("GenomicRanges")
a



#needed to install rtracklayer

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")  # Install BiocManager if it's not already installed
}
BiocManager::install("rtracklayer")



library(GenomicRanges)
library(rtracklayer) # for reading in bed file
library(here)
library(usethis)
library(devtools)

TAZ_peaks <- import(here("~/Desktop/Computational_Genomics/Midterm Presentation/TAZ_peak/TAZ_peaks.narrowPeak"))
YAP_peaks <- import(here("~/Desktop/Computational_Genomics/Midterm Presentation/YAP_peak/YAP_peaks.narrowPeak"))

TAZ_peaks
YAP_peaks

seqlevels(TAZ_peaks)
seqlevels(YAP_peaks)

common_levels <- intersect(seqlevels(TAZ_peaks), seqlevels(YAP_peaks))

TAZ_peaks <- keepSeqlevels(TAZ_peaks, common_levels, pruning.mode = "coarse")
YAP_peaks <- keepSeqlevels(YAP_peaks, common_levels, pruning.mode = "coarse")

TAZ_overlap_YAP_peaks<- subsetByOverlaps(TAZ_peaks, YAP_peaks)
length(TAZ_overlap_YAP_peaks)

YAP_overlap_TAZ_peaks<- subsetByOverlaps(YAP_peaks, TAZ_peaks)
length(YAP_overlap_TAZ_peaks)



#download package from github
#install_github("js229/Vennerable")
#giving PAT error


#install these packages to use Vennerable
BiocManager::install(c("graph", "RBGL"))

#used git clone https://github.com/js229/Vennerable.git
#in terminal to copy to local
#cd into Vennerable and pwd
#install Vennerable in R:
remotes::install_local("/Users/veronicagosnell/Vennerable")

library(Vennerable)

n_YAP <- length(YAP_peaks)  # Total peaks 
n_TAZ <- length(TAZ_peaks)  # Total peaks 

n_overlap <- length(YAP_overlap_TAZ_peaks)

venn_data <- Venn(SetNames = c("YAP", "TAZ"),
                  Weight = c(
                    "10" = n_YAP, # Unique to A
                    "01" = n_TAZ, # Unique to B
                    "11" = n_overlap         # Intersection
                  ))
n_TAZ
# Plot the Venn diagram
plot(venn_data)

TEAD4_peak<- import(here("~/Desktop/Computational_Genomics/Midterm Presentation/TEAD4_peak/TEAD4_peaks.narrowPeak"))

YAP_overlap_TAZ_peaks_overlap_TEAD4<- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peak)

n_YAP_TAZ <- length(YAP_overlap_TAZ_peaks)  # Total peaks 
n_TEAD4 <- length(TEAD4_peak)  # Total peaks 
n_overlap2<- length(YAP_overlap_TAZ_peaks_overlap_TEAD4)

venn_data2 <- Venn(SetNames = c("YAP/TAZ", "TEAD4"),
                   Weight = c(
                     "10" = n_YAP_TAZ, # Unique to A
                     "01" = n_TEAD4, # Unique to B
                     "11" = n_overlap2        # Intersection
                   ))

# Plot the Venn diagram
plot(venn_data2)
