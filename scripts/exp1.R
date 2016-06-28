source("../scripts/multi_index.R")
source("../external/sdsl-lite/benchmark/basic_functions.R")
library(ggplot2) 

options(scipen=999)

data <- data_frame_from_key_value_pairs("../build/results/exp1.result.txt")

data[["hash_file"]] <- pretty_input(data[["hash_file"]])


# Check that for each (k, hash_file, hashes, queries)-combination the
# checksum is the same

dd <- data[c("index","k","hash_file","hashes","queries","check_unique_matches_full")]

if ( sum(duplicated(unique(dd[-1])[-5])) > 0 ){
    pdf("exp1.pdf")
    plot(c(),c(),xlim=c(0,1),ylim=c(0,0.1),axes=F,ylab=NA,xlab=NA,main="checksum error in exp1")    
    dev.off()
    q()
}

d <- data[c("hash_file","index","qry_file","b","k","time_per_full_query_in_us","time_per_search_query_in_us")] 

colnames(d) <- c("Input","Index","Queries","b","e","time","stime")

d[["etime"]] <- d[["time"]]-d[["stime"]]
d[["Candidates"]] <- data[["candidates_per_query"]]

d[["Queries"]] <- basename(as.character(data[["qry_file"]]))
d[["Queries"]] <- gsub("\\.query$","",d[["Queries"]])
d[["Queries"]] <- gsub("(.*\\.)(.*)","\\2",d[["Queries"]])
d <- cbind(d,"size_per_key" = data[["index_size_in_bytes"]]/data[["hashes"]])

d[["Index"]] <- gsub("mi_(.*)_red$","mi_red_\\1",d[["Index"]])
d[["Method"]] <- gsub(".*_(.*)$","\\1",d[["Index"]])
d[["Index"]] <- gsub("(.*)_(.*)$","\\1",d[["Index"]])

dd <- d[c("size_per_key","Index","Queries","e","Method","time","stime","etime","Input")]
ddd <- dd

ddd[["Index"]] <- gsub("_","\\\\_",ddd[["Index"]])
ddd[["Method"]] <- gsub("bv","succinct",ddd[["Method"]])
ddd[["Method"]] <- gsub("bs","bin. search",ddd[["Method"]])

d4 <- ddd[ddd$e %in% c(3,4),]
d4 <- subset(d4,d4$Method %in% c("bin. search","succinct"))
d4 <- d4[grep("red",d4$Index),]

d4$stimepercent <- d4[["stime"]]/d4[["time"]]*100

d4c <- rbind(
cbind(d4,"part"=rep("searching",nrow(d4)),"parttime"=d4$stime) ,
cbind(d4,"part"=rep("checking",nrow(d4)),"parttime"=d4$etime) 
)

d4c$k <- gsub("(.*)","$k=\\1$",as.character(d4c$e))
d4c$check_percentage <- sprintf("(%.0f \\%%)", d4c$stime/d4c$time*100)
d4c[d4c$part=="checking",]$check_percentage <- ""

plot <- qplot(Method,
             data=d4c,
             fill=part,
             weight=parttime,
             xlab="Search method",
             ylab="Average query time [$\\mu s$]"
             ) +
        facet_grid(k ~ Input) +
        scale_fill_discrete(guide = guide_legend(title = "Query phase", reverse = TRUE)) +
        geom_text(aes(y=parttime + 0.03*max(time), ymax=max(parttime), label=check_percentage),
                  data=d4c, size=4,
                  color= scales::hue_pal()(3)[1],
                  vjust=0
                 )
     

pdf("exp1.pdf")
print(plot)
dev.off()
