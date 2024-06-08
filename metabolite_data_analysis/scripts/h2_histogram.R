
library(data.table)
library(ggplot2)

file <- "data/tab_h2.txt"

h2_tmp <- fread(file)

h2 <- h2_tmp[h2_tmp$Type == "single_GRM",]

mean_h2 <- mean(h2$Vg_Vp)

pdf("results/metab_h2_histogram.pdf")
ggplot(h2,aes(x = Vg_Vp)) +
        geom_histogram(binwidth=0.01) +
        xlab("Heritability Estimate") +
        ylab("Number of Metabolites") +
        geom_vline(xintercept = mean_h2,color = "red",alpha=0.5) +
        annotate("text", label = paste0("Mean = ",round(mean_h2,digits=2)), x = 0.4, y = 70, alpha = 0.5,size = 8, colour = "red") +
        theme_bw() +
        theme(text = element_text(size=10,face="bold"),aspect.ratio=1)
dev.off()

