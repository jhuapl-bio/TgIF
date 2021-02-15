library(ggplot2)

## string vector containing entries from command line (or bash wrapper)
args = commandArgs(trailingOnly=TRUE)
input<-args[1];	# position.data file

#input <- '/home/playera1/APL/Shared_Documents/Documents/DARPA-LiSTENS/fy20/transgene_stuff/tgif_algorithm_dev_3cg2/tgif-20200531_FLOWCELLID_3cg2_shear_fastq_pass.fastq/plots/Chr08_33451539_33451905.data'

d <- read.table(input, header=F, sep='\t')



filename <- basename(input)
outdir <- dirname(input)



df <- d
pos <- subset(df, V3>0)
pos$length <- 1000000000-(pos$V4-pos$V2)
neg <- subset(df, V3<0)
neg$length <- (neg$V4-neg$V2)-1000000000
df <- rbind(pos,neg)
# order legend
df$V1 <- factor(df$V1, levels = c("(+)","(-)"))
# x-axis scale
xmin <- round(min(df$V2),digits=-3)
xmax <- round(max(df$V4),digits=-3)
b <- (xmax-xmin)/10

t <- theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y = element_blank())
labels <- labs(x="genome position (1-based)", y="", color="strand")
p <- ggplot(df) + geom_segment(aes(x=V2, y=reorder(V3, length), xend=V4, yend=reorder(V3, length), color=V1), alpha=0.5, size=0.75) + scale_x_continuous(breaks=seq(xmin, xmax, b)) + scale_color_manual(values = c("red","blue")) + labels + guides(size = FALSE, alpha = FALSE) + t

# adjust height base on total reads aligned (nrows)
h=1+(nrow(d)/25)
outfile=paste(outdir,"/",filename,"_gapmap.png",sep="")
ggsave(outfile, p, width=4, height=h, units="in", device="png", dpi="retina", type="cairo")





