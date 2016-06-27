require("ggplot2")

if ( !exists( "tikzDeviceLoaded" ) ){  
	require(tikzDevice) #if not installed call install.packages("tikzDevice", repos="http://R-Forge.R-project.org")
	options(tikzLatexPackages = c(getOption('tikzLatexPackages'),
        	paste("\\input{",getwd(),"/../paper/defs.tex}",sep="")))
	tikzDeviceLoaded = T
}

args <- commandArgs(trailingOnly = TRUE)
info_file <- args[1]

files <- list.files("../data/", ".distr")
files <- gsub(".data.distr$","",files)

df <- data.frame("x" = seq(0,64), "y" = choose(64, seq(0,64))/2**64, "Input" = "uniform" )

for(f in files){
    d <- read.csv(paste("../data/",f,".data.distr",sep=""),header=F)    
    colnames(d) <- c("x","y")
    d["y"] <- d["y"]/sum(d["y"])
    d["Input"] <- rep(f,65)
    df <- rbind(df, d)
}

df[["Input"]] <- gsub("Clueweb09-Full.SimHash","\\\\cluesim",df[["Input"]])
df[["Input"]] <- gsub("Clueweb09-Full.OddSketch","\\\\clueodd",df[["Input"]])
df[["Input"]] <- gsub("lsh_sift_64.hash","\\\\siftlsh",df[["Input"]])
df[["Input"]] <- gsub("mlh_sift_64.hash","\\\\siftmlh",df[["Input"]])

tikz('input_info.tex', standAlone = FALSE, width=3.5, height=3.2)

plot <- qplot(x,y,data=df,
     color=Input,
#     shape=Input,
     xlab="Hamming weight",
     ylab="Fraction of keys"
#     ,log="y"
     ) + theme_bw(base_size = 12, base_family = "") + 
      theme(legend.position="top") + geom_line() +
      guides(color=guide_legend(ncol=3)) +
      scale_color_discrete("") +
      scale_x_continuous(breaks=c(0,16,32,48,64))

print(plot)

dev.off()
