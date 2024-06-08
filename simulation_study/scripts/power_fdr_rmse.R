#Power and FDR tables

library(qvalue)
library(ggplot2)
library(INTACT)
library(dplyr)
library(tidyr)
library(data.table)

fdr_rst <- function (posterior, alpha = 0.05)
{
    gene_num <- seq(1, length(posterior))
    lfdr_sort = sort(1 - posterior)
    FDR = cumsum(lfdr_sort)/(1:length(lfdr_sort))
    thresh = 1 - lfdr_sort[max(which(FDR <= alpha))]
    rej_gene = as.numeric(gene_num[which(posterior > thresh)])
    out_tmp <- rep(FALSE, length(posterior))
    out_tmp[rej_gene] <- TRUE
    out <- data.frame(posterior = posterior, sig = out_tmp)
    return(out)
}

d <- fread("sim_data/rst_summary.txt",h=T)


d$sim_gene = paste0(d$sim_num,"_",d$gene)

cgene = d$sim_gene[which(d$true_effect != 0)]

d$intact_ptwas_50 = intact(z_vec = d$ptwas_z,GLCP_vec = d$GRCP,t=0.5)

power_fdr <- function(stat_name,type,stat_effect_name=NULL){

        posteriors <- d[[paste0(stat_name)]]

        if(type == "posterior"){

        lfdr <- sort(1-posteriors)

        FDR <- cumsum(lfdr)/(1:length(lfdr))

        thresh = 1-lfdr[max(which(FDR<=0.05))]

        rej = d$sim_gene[posteriors>thresh]

        }

        if(type == "z"){

        p = 2*pnorm(-abs(d[[paste0(stat_name)]]))

        qrst = qvalue(p,fdr.level=0.05)

        rej=d$sim_gene[which(qrst$sig)]

        }

        fp_set = setdiff(rej, cgene)

        tp_set = intersect(rej, cgene)

        FDR = round(length(fp_set)/length(rej),digits=3)

        power = round(length(tp_set)/length(cgene),digits=3)

        out = data.frame("method" = stat_name,
                          "n_rej" = length(rej),
                          "Power" = power,
                          "FDR" = FDR)
        if(!is.null(stat_effect_name)){

        rmse = round(sqrt(mean((d[[paste0(stat_effect_name)]] - d$true_effect)^2)),digits=3)

        out$RMSE = rmse

        }

        return(out)
}

#Marginal methods

out_marginal <- rbind.data.frame(power_fdr("ptwas_z","z"),
                        power_fdr("elastic_net_z","z"),
                        power_fdr("smr_z","z"),
                        power_fdr("GLCP","posterior"),
                        power_fdr("GRCP","posterior"),
                        power_fdr("intact_ptwas_50","posterior"),
			power_fdr("susie_pip_pgene_grcp_dap_25_dap","posterior"))

out_marginal$method <- c("PTWAS","PrediXcan","SMR","GLCP","GRCP","INTACT","ISuSiE")

print(out_marginal)

#Joint methods

out_joint <- rbind.data.frame(power_fdr("susie_pip_multi_gene_twas_dap_25","posterior","susie_effect_multi_gene_twas_dap_25"),                       
			      power_fdr("susie_pip_default_dap_25","posterior","susie_effect_default_dap_25"),                       
			      power_fdr("susie_pip_focus_dap_25","posterior","susie_effect_focus_dap_25"),                        
			      power_fdr("susie_pip_ctwas_priors_ptwas","posterior","susie_effect_ctwas_priors_ptwas"),
			      power_fdr("susie_pip_pgene_grcp_dap_25_dap","posterior","susie_effect_pgene_grcp_dap_25_dap"), 
			      power_fdr("susie_pip_true_priors_censored_wts_ptwas","posterior","susie_effect_true_priors_censored_wts_ptwas"))



out_joint$method <- c("Multi-gene TWAS","SuSiE Default","FOCUS","cTWAS","ISuSiE","Oracle Priors")

print(out_joint)


#Vary ISuSiE TWAS prediction models

out_pred_model <- rbind.data.frame(power_fdr("susie_pip_pgene_grcp_dap_25_dap","posterior","susie_effect_pgene_grcp_dap_25_dap"),
				   power_fdr("susie_pip_pgene_grcp_dap_25_elastic_net","posterior","susie_effect_pgene_grcp_dap_25_elastic_net"),
				   power_fdr("susie_pip_pgene_grcp_dap_25_smr","posterior","susie_effect_pgene_grcp_dap_25_smr"))

out_pred_model$method <- c("PTWAS","PrediXcan","SMR")

print(out_pred_model)


#Vary CPIP threshold

out_cpip <- rbind.data.frame(power_fdr("susie_pip_pgene_grcp_dap_25_dap","posterior","susie_effect_pgene_grcp_dap_25_dap"),
			     power_fdr("susie_pip_pgene_grcp_dap_50_dap","posterior","susie_effect_pgene_grcp_dap_50_dap"),
			     power_fdr("susie_pip_pgene_grcp_dap_75_dap","posterior","susie_effect_pgene_grcp_dap_75_dap"),
			     power_fdr("susie_pip_pgene_grcp_dap_95_dap","posterior","susie_effect_pgene_grcp_dap_95_dap"))

colnames(out_cpip)[1] <- "CPIP threshold"

out_cpip[["CPIP threshold"]] <- c(0.25,0.50,0.75,0.95)

print(out_cpip)


#Use SuSiE credible sets for instrument filtering

out_susie_cs <- rbind.data.frame(power_fdr("susie_pip_pgene_grcp_25_susie","posterior","susie_effect_pgene_grcp_25_susie"),
                        power_fdr("susie_pip_pgene_grcp_50_susie","posterior","susie_effect_pgene_grcp_50_susie"),
                        power_fdr("susie_pip_pgene_grcp_75_susie","posterior","susie_effect_pgene_grcp_75_susie"),
                        power_fdr("susie_pip_pgene_grcp_95_susie","posterior","susie_effect_pgene_grcp_95_susie"))

colnames(out_susie_cs)[1] <- "susie credible set coverage"

out_susie_cs[["susie credible set coverage"]] <- c("0.25","0.50","0.75","0.95")


#Power & FDR vs CPIP threshold plot

out_thresh <- rbind.data.frame(power_fdr("susie_pip_pgene_grcp_pi1_ptwas","posterior","susie_effect_pgene_grcp_pi1_ptwas"),
                        power_fdr("susie_pip_pgene_grcp_dap_25_dap","posterior","susie_effect_pgene_grcp_dap_25_dap"),
                        power_fdr("susie_pip_pgene_grcp_dap_50_dap","posterior","susie_effect_pgene_grcp_dap_50_dap"),
                        power_fdr("susie_pip_pgene_grcp_dap_75_dap","posterior","susie_effect_pgene_grcp_dap_75_dap"),
                        power_fdr("susie_pip_pgene_grcp_dap_95_dap","posterior","susie_effect_pgene_grcp_dap_95_dap"))

out_thresh$thresh <- c(0,0.25,0.5,0.75,0.95)

pdf("sim_rst/power_fdr_cpip_threshold.pdf")
out_thresh %>%
        pivot_longer(cols = Power:FDR,
                     names_to = "type",
                     values_to = "stat") %>%
        ggplot(aes(x = thresh,y=stat,col=type)) +
        geom_point() +
        xlab("CPIP Threshold") +
        ylab("Statistic") +
        geom_line() +
        theme_bw() +
        theme(text = element_text(size = 10,face="bold"),legend.title=element_blank(),aspect.ratio=1)
dev.off()






