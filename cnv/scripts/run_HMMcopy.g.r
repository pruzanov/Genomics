# Below commands are for creating segmentation files and producing some QC images
# outputBasename used to construct names for files produced by this script

library(HMMcopy)

cmd_args=commandArgs(trailingOnly = TRUE)
# Arguments should be passed as: normal.wig, tumor.wig, refGC.wig, ref_mappable.wig, outputBasename
normalReads<-cmd_args[1]
gcContent<-cmd_args[2]
refMappable<-cmd_args[3]
outputBasename<-cmd_args[4]

norm_uncorrected_reads <- wigsToRangedData(normalReads, gcContent, refMappable)
norm_corrected_copy <- correctReadcount(norm_uncorrected_reads)

# Should take no longer than a few minutes on a human genome.
# The correctReadcount requires at least about 1000 bins to work properly.
# 2. Segmentation

# Below commands in R
param <- HMMsegment(norm_corrected_copy, getparam = TRUE) # retrieve converged parameters via EM
param$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2)
param$m <- param$mu
segmented_copy <- HMMsegment(norm_corrected_copy, param) # perform segmentation via Viterbi

# 3. Export
# Export to SEG format for CNAseq segmentation
segFile<-paste(outputBasename, "seg", sep = ".")
tsvFile<-paste(outputBasename, "tsv", sep = ".")

rangedDataToSeg(norm_corrected_copy, file = segFile)
write.table(segmented_copy$segs, file = tsvFile, quote = FALSE, sep = "\t")

options(bitmapType="cairo")
# 4. Visualization - produce some images with hard-coded dimensions

# Bias plots:
print("Producing CG Bias plot...")
png(filename = paste(outputBasename,"bias_plot","png", sep="."),width = 1200, height = 580, units="px", pointsize=15, bg="white")
plotBias(norm_corrected_copy)  # May be one plot per comparison  1200x580
dev.off()

chroms<-unique(segmented_copy$segs$chr)

# Segmentation plots:
# need to do it one plot per chromosome 1200x450
print("Producing Segmentation plots...")
par(mfrow = c(1, 1))
for (c in 1:length(chroms)) {
 if (!grepl("_",chroms[c]) && !grepl("M",chroms[c])) {
	 png(filename = paste(outputBasename,"s_plot", chroms[c], "png", sep="."),width = 1200, height = 450, units="px", pointsize=15, bg="white")
	 plotSegments(norm_corrected_copy, segmented_copy, pch = ".", ylab = "Germline CNV results", xlab = "Chromosome Position",chr = chroms[c], main = paste("Segmentation for Chromosome",chroms[c], sep=" "))
	 cols <- stateCols()
	 legend("topleft", c("HOMD", "HETD", "NEUT", "GAIN", "AMPL", "HLAMP"), fill = cols, horiz = TRUE, bty = "n", cex = 0.9)
	 dev.off()
 } else {
	 print(paste("Chromosome",c,"cannot be plotted with  plotSegments",sep = " "))
 }
}


# Correction plots:
# need to do it one plot per chromosome 1200x680

print("Producing Correction plots...")


for (c in 1:length(chroms)) {
 
 if (!grepl("_",chroms[c]) && !grepl("M",chroms[c])) {
  print(paste("Cleaning data for",chroms[c],sep=" "))
  
  slice<-norm_corrected_copy[paste(chroms[c])]
  slice$sd<-slice$reads/median(slice$reads, na.rm = TRUE)
  slice.clean<-slice[!is.na(slice$sd),]
  
  print(paste("Plotting for chromosome",c,sep=" "))
  png(filename = paste(outputBasename,"c_plot",chroms[c],"png", sep="."),width = 1200, height = 680, units="px", pointsize=15, bg="white")
  plotCorrection(slice.clean, pch = ".", chr= chroms[c], na.rm = TRUE)
  dev.off()
 } else {
  print(paste("Chromosome",c,"cannot be plotted with  plotCorrection",sep = " "))
 }
}

