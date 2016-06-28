library(ggplot2)
source("../scripts/multi_index.R")
source("../external/sdsl-lite/benchmark/basic_functions.R")

options(scipen=999)

data <- data_frame_from_key_value_pairs("../build/results/exp0.result.txt")
data[["hash_file"]] <- pretty_input(data[["hash_file"]])
data[["time"]] <- data[["time_per_full_query_in_us"]] #/ ( data[["check_unique_matches_full"]] +1 )
data[["query"]] <- data[["qry_file"]]
d <- data

# Check that for each (k, hash_file, hashes, queries)-combination the
# checksum is the same

dd <- d[c("index","k","hash_file","hashes","queries","check_unique_matches_full")]

if ( sum(duplicated(unique(dd[-1])[-5])) > 0 ){
    pdf("exp0.pdf")
    plot(c(),c(),xlim=c(0,1),ylim=c(0,0.1),axes=F,ylab=NA,xlab=NA,main="checksum error in exp0")    
    dev.off()
    q()
}

# Generate plot
plot <- ggplot(d, aes(x=factor(k), y=time, fill=index),
              ylab="Average query time [$\\mu s$]",
              xlab="$k$"
              ) +
         geom_bar(position='dodge',stat='identity') +
        scale_y_continuous(name="Average query time [$\\mu s$]", limits=c(1,2*max(d$time)), trans='log10') +
        scale_x_discrete("k", factor(d$k)) +
        facet_wrap(~hash_file) +
        scale_fill_discrete(guide = guide_legend(title="Index")  ) +
        geom_text(aes(y=2*max(time), ymax=2*max(time), label=sprintf("(%.1f)",index_size_in_bytes/(8*hashes)), color=index),
                  position = position_dodge(width=1),
                  data=d, size=3  ) +
        scale_color_discrete(guide = guide_legend(title="Index")  ) +
        geom_text(aes(y=time/1.5, ymax=max(time), label=sprintf("%d",round(time))),
                  position = position_dodge(width=1),
                  data=d, size=3,#,angle=-90
                  color="white"
        )
pdf("exp0.pdf")
print(plot)
dev.off()



