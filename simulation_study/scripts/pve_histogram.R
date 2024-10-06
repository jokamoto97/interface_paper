#Make histogram to visualize distribution of PVE for GWAS simulation

library(ggplot2)
library(dplyr)

pve <- read.table("sim_data/gwas_pve.txt")

mean_pve <- mean(pve$V1)

pdf("sim_rst/gwas_pve.pdf")
ggplot(pve,aes(x = V1)) +
        geom_histogram(binwidth=0.01) +
        xlab("PVE") +
        ylab("Number of Simulated GWAS Traits") +
        geom_vline(xintercept = mean_pve,color = "red",alpha=0.5) +
        annotate("text", label = paste0("Mean PVE = ",round(mean_pve,digits=2)), x = 0.3, y = 40, alpha = 0.5,size = 8, colour = "red") +
        theme_bw() +
        theme(text = element_text(size=10,face="bold"),aspect.ratio=1)
dev.off()

#Make histogram to visualize distribution of PVE for eQTL simulation

pve <- read.table("sim_data/eqtl.truth") %>%
        dplyr::select(V1,V4) %>%
        unique()

mean_pve <- mean(pve$V4)

#600 region datasets
pve <- do.call("rbind", replicate(600, pve, simplify = FALSE))

pdf("sim_rst/expr_pve.pdf")
pve %>%
        ggplot(aes(x = V4)) +
        geom_histogram(binwidth=0.05) +
        xlab("PVE") +
        ylab("Number of Simulated Genes") +
        geom_vline(xintercept = mean_pve,color = "red",alpha=0.5) +
        annotate("text", label = paste0("Mean Expression PVE = ",round(mean_pve,digits=2)), x = 0.3, y = 5000, alpha = 0.5,size = 6, colour = "red") +
        theme_bw() +
        theme(text = element_text(size=10,face="bold"),aspect.ratio=1)
dev.off()

