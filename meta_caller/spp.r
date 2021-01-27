library(spp)
library(snow)

filenames <- commandArgs(trailingOnly=TRUE)
file.data <- filenames[1]
file.input <- filenames[2]
prefix <- filenames[3]
cluster <- makeCluster(as.numeric(filenames[4]))

input.data <- read.bam.tags(file.input, read.tag.names=TRUE)
chip.data <- read.bam.tags(file.data, read.tag.names=TRUE)

setwd("spp/")

path <- getwd()

#get binding info from cross-correlation profile
binding.characteristics <- get.binding.characteristics(chip.data,srange=c(50,500),bin=5, cluster=cluster, accept.all.tags=FALSE)

# Print out binding peak separation distance
print(paste("binding peak separation distance =",binding.characteristics$peak$x))

# Plot cross-correlation profile
pdf(file=paste(path,"/", prefix, ".crosscorrelation.pdf", sep=""),width=5,height=5)
par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation")
abline(v=binding.characteristics$peak$x,lty=2,col=2)
dev.off()

# select informative tags based on the binding characteristics
chip.data <- select.informative.tags(chip.data, binding.characteristics)
input.data <- select.informative.tags(input.data, binding.characteristics)

#Subtract background
# restrict or remove singular positions with very high tag counts
chip.data <- remove.local.tag.anomalies(chip.data)
input.data <- remove.local.tag.anomalies(input.data)

#Determine binding positions
# binding detection parameters
# desired FDR (1%). Alternatively, an E-value can be supplied to the method calls below instead of the fdr parameter
fdr <- 0.9

# the binding.characteristics contains the optimized half-size for binding detection window
detection.window.halfsize <- binding.characteristics$whs

# determine binding positions using wtd method
bp <- find.binding.positions(signal.data=chip.data,control.data=input.data,fdr=fdr,whs=detection.window.halfsize,cluster=cluster)

print(paste("detected",sum(unlist(lapply(bp$npl,function(d) length(d$x)))),"peaks"))


# output detected binding positions
output.binding.results(bp,paste(path, "/", prefix,".binding.positions.txt", sep=""))
write.narrowpeak.binding(bp,paste(path, "/", prefix,".narrowPeak", sep=""))
