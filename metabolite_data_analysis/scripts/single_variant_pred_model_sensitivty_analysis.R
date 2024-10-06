#Plot effect size heterogeneity when using various single-SNP prediction models

library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)


get_i2 <- function(betas,sds){

        wi <- 1/(sds^2)

        mean_beta <- sum(wi*betas)/sum(wi)

        Q <- sum((betas - mean_beta)^2 / sds^2)

        t <- length(betas)

        i2 <- max(0,(Q - t + 1)/Q)

        return(i2)
}


rst <- fread("data/smr_diff_set.txt",h=TRUE,data.table=F)

if(!dir.exists("results/intact_interface_diff_set_effect_size_forest_plots_smr/")){dir.create("results/intact_interface_diff_set_effect_size_forest_plots_smr/")}

make_forest_plot <- function(cid){

        rst_tmp <- rst[rst$cid == cid,]

        gene <- unique(rst_tmp$gene)

        rst_tmp$Top_pQTL_CS_ID <- factor(rst_tmp$Top_pQTL_CS_ID,levels = rev(rst_tmp$Top_pQTL_CS_ID))

        rst_tmp$upper <- rst_tmp$Beta + 1.96*rst_tmp$SE

        rst_tmp$lower <- rst_tmp$Beta - 1.96*rst_tmp$SE

        I2 <- round(get_i2(rst_tmp$Beta,rst_tmp$SE),digits = 3)

        pdf(paste0("results/intact_interface_diff_set_effect_size_forest_plots_smr/",cid,"_",gene,"_forest_plot.pdf"))
        print(
              ggplot(rst_tmp,aes(x = Top_pQTL_CS_ID,y=Beta,ymin=lower,ymax=upper)) +
                      geom_pointrange() +
                      geom_hline(yintercept = 0,lty = 2) +
                      coord_flip() +
                      xlab("Protein Prediction Model Credible Set ID") +
                      ylab("Gene-to-trait Effect (95% CI)") +
                      annotate(geom = 'text',
                               label =paste0("I\u00b2  = ", I2),
                               x = Inf, y = -Inf, hjust = -0.2, vjust = 2) +
                      theme_bw() +
                      theme(text = element_text(size=15,face="bold"),aspect.ratio=1)
        )
        dev.off()
}


cids <- unique(rst$cid)

lapply(cids,make_forest_plot)



