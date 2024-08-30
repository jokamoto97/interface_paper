#Gene-to-trait effect estimate vs true effect scatterplots

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


plot_effects <- function(stat_names,effect_names,facet_labels,filename,sd_names=NULL){

        n_stats = length(stat_names)

        for (i in 1:n_stats){

        d[[paste0(stat_names[i],"_sig")]] <- fdr_rst(d[[stat_names[i]]])$sig

        }

        if(!is.null(sd_names)){

                for (i in 1:n_stats){
                d[[paste0(effect_names[i],"_upper")]] = d[[effect_names[i]]] + 1.96*d[[sd_names[i]]]
                d[[paste0(effect_names[i],"_lower")]] = d[[effect_names[i]]] - 1.96*d[[sd_names[i]]]
                }
        }
        plotdf <- d %>%
                pivot_longer(cols = all_of(effect_names),
                             names_to="Type",
                             values_to = "susie_effect") %>%
                mutate(Type = case_when(Type == effect_names[1] ~ facet_labels[1],
                                        Type == effect_names[2] ~ facet_labels[2],
                                        Type == effect_names[3] ~ facet_labels[3]))
        plotdf$susie_sig = FALSE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[1] & plotdf[[paste0(stat_names[1],"_sig")]] == TRUE] = TRUE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[2] & plotdf[[paste0(stat_names[2],"_sig")]] == TRUE] = TRUE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[3] & plotdf[[paste0(stat_names[3],"_sig")]] == TRUE] = TRUE
#        plotdf <- plotdf %>% filter(susie_sig == TRUE)

        if(!is.null(sd_names)){
                plotdf$upper_lim = NA
                plotdf$lower_lim = NA
                for (i in 1:n_stats){
                        plotdf$upper_lim[plotdf$Type == facet_labels[i]] = plotdf[[paste0(effect_names[i],"_upper")]][plotdf$Type == facet_labels[i]]
                        plotdf$lower_lim[plotdf$Type == facet_labels[i]] = plotdf[[paste0(effect_names[i],"_lower")]][plotdf$Type == facet_labels[i]]
                }
                #only show for outliers
                plotdf$upper_lim[plotdf$upper_lim - plotdf$susie_effect < 5] = 0
                plotdf$lower_lim[plotdf$susie_effect - plotdf$lower_lim < 5] = 0
        }
        if(is.null(sd_names)){
        pdf(paste0("sim_rst/",filename))
        print(
        plotdf %>%
		mutate(across(Type, ~factor(., levels=rev(facet_labels)))) %>%
                ggplot(aes(x = true_effect, y = susie_effect)) +
                geom_point() +
                geom_abline(intercept = 0, slope = 1, color = "red") +
                #coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
                facet_wrap(~Type,scales = "free") +
                xlab("True gene-to-trait effect size") +
                ylab("SuSiE regression coefficient estimate") +
                theme_bw() +
                theme(text = element_text(size = 10, face = "bold"),aspect.ratio=1)
        )
        dev.off()
        }
        if(!is.null(sd_names)){
        pdf(paste0("sim_rst/",filename))
        print(
        plotdf %>%
                mutate(across(Type, ~factor(., levels=c(facet_labels[2],facet_labels[3],facet_labels[1])))) %>%
                ggplot(aes(x = true_effect, y = susie_effect)) +
                geom_point() +
                geom_errorbar(aes(ymin=lower_lim, ymax=upper_lim)) +
                geom_abline(intercept = 0, slope = 1, color = "red") +
                #coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
                facet_wrap(~Type) +
                xlab("True gene-to-trait effect size") +
                ylab("SuSiE regression coefficient estimate") +
                theme_bw() +
                theme(text = element_text(size = 10, face = "bold"),aspect.ratio=1)
        )
        dev.off()
        }
}


plot_effects2 <- function(stat_names,effect_names,facet_labels,filename,sd_names=NULL){

        n_stats = length(stat_names)

        for (i in 1:n_stats){

        d[[paste0(stat_names[i],"_sig")]] <- fdr_rst(d[[stat_names[i]]])$sig

        }

        if(!is.null(sd_names)){

                for (i in 1:n_stats){
                d[[paste0(effect_names[i],"_upper")]] = d[[effect_names[i]]] + 1.96*d[[sd_names[i]]]
                d[[paste0(effect_names[i],"_lower")]] = d[[effect_names[i]]] - 1.96*d[[sd_names[i]]]
                }
        }
        plotdf <- d %>%
                pivot_longer(cols = all_of(effect_names),
                             names_to="Type",
                             values_to = "susie_effect") %>%
                mutate(Type = case_when(Type == effect_names[1] ~ facet_labels[1],
                                        Type == effect_names[2] ~ facet_labels[2],
                                        Type == effect_names[3] ~ facet_labels[3],
                                        Type == effect_names[4] ~ facet_labels[4]))
        plotdf$susie_sig = FALSE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[1] & plotdf[[paste0(stat_names[1],"_sig")]] == TRUE] = TRUE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[2] & plotdf[[paste0(stat_names[2],"_sig")]] == TRUE] = TRUE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[3] & plotdf[[paste0(stat_names[3],"_sig")]] == TRUE] = TRUE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[4] & plotdf[[paste0(stat_names[4],"_sig")]] == TRUE] = TRUE
 #       plotdf <- plotdf %>% filter(susie_sig == TRUE)

        if(!is.null(sd_names)){
                plotdf$upper_lim = NA
                plotdf$lower_lim = NA
                for (i in 1:n_stats){
                        plotdf$upper_lim[plotdf$Type == facet_labels[i]] = plotdf[[paste0(effect_names[i],"_upper")]][plotdf$Type == facet_labels[i]]
                        plotdf$lower_lim[plotdf$Type == facet_labels[i]] = plotdf[[paste0(effect_names[i],"_lower")]][plotdf$Type == facet_labels[i]]
                }
                #only show for outliers
                plotdf$upper_lim[plotdf$upper_lim - plotdf$susie_effect < 5] = 0
                plotdf$lower_lim[plotdf$susie_effect - plotdf$lower_lim < 5] = 0
        }
        if(is.null(sd_names)){
        pdf(paste0("sim_rst/",filename))
        print(
        plotdf %>%
                ggplot(aes(x = true_effect, y = susie_effect)) +
                geom_point() +
                geom_abline(intercept = 0, slope = 1, color = "red") +
                #coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
                facet_wrap(~Type,scales = "free") +
                xlab("True gene-to-trait effect size") +
                ylab("SuSiE regression coefficient estimate") +
                theme_bw() +
                theme(text = element_text(size = 10, face = "bold"),aspect.ratio=1)
        )
        dev.off()
        }
        if(!is.null(sd_names)){
        pdf(paste0("sim_rst/",filename))
        print(
       plotdf %>%
                ggplot(aes(x = true_effect, y = susie_effect)) +
                geom_point() +
                geom_errorbar(aes(ymin=lower_lim, ymax=upper_lim)) +
                geom_abline(intercept = 0, slope = 1, color = "red") +
                #coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
                facet_wrap(~Type,scales = "free") +
                xlab("True gene-to-trait effect size") +
                ylab("SuSiE regression coefficient estimate") +
                theme_bw() +
                theme(text = element_text(size = 10, face = "bold"),aspect.ratio=1)
        )
        dev.off()
        }
}


plot_effects3 <- function(stat_names,effect_names,facet_labels,filename,sd_names=NULL){

        n_stats = length(stat_names)

        for (i in 1:n_stats){

        d[[paste0(stat_names[i],"_sig")]] <- fdr_rst(d[[stat_names[i]]])$sig

        }

        if(!is.null(sd_names)){

                for (i in 1:n_stats){
                d[[paste0(effect_names[i],"_upper")]] = d[[effect_names[i]]] + 1.96*d[[sd_names[i]]]
                d[[paste0(effect_names[i],"_lower")]] = d[[effect_names[i]]] - 1.96*d[[sd_names[i]]]
                }
        }
        plotdf <- d %>%
                pivot_longer(cols = all_of(effect_names),
                             names_to="Type",
                             values_to = "susie_effect") %>%
                mutate(Type = case_when(Type == effect_names[1] ~ facet_labels[1],
                                        Type == effect_names[2] ~ facet_labels[2],
                                        Type == effect_names[3] ~ facet_labels[3],
                                        Type == effect_names[4] ~ facet_labels[4]#,
         #                               Type == effect_names[5] ~ facet_labels[5]
			     )
                )
        plotdf$susie_sig = FALSE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[1] & plotdf[[paste0(stat_names[1],"_sig")]] == TRUE] = TRUE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[2] & plotdf[[paste0(stat_names[2],"_sig")]] == TRUE] = TRUE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[3] & plotdf[[paste0(stat_names[3],"_sig")]] == TRUE] = TRUE
        plotdf$susie_sig[plotdf[["Type"]] == facet_labels[4] & plotdf[[paste0(stat_names[4],"_sig")]] == TRUE] = TRUE
        #plotdf$susie_sig[plotdf[["Type"]] == facet_labels[5] & plotdf[[paste0(stat_names[5],"_sig")]] == TRUE] = TRUE
  #      plotdf <- plotdf %>% filter(susie_sig == TRUE)

        if(!is.null(sd_names)){
                plotdf$upper_lim = NA
                plotdf$lower_lim = NA
                for (i in 1:n_stats){
                        plotdf$upper_lim[plotdf$Type == facet_labels[i]] = plotdf[[paste0(effect_names[i],"_upper")]][plotdf$Type == facet_labels[i]]
                        plotdf$lower_lim[plotdf$Type == facet_labels[i]] = plotdf[[paste0(effect_names[i],"_lower")]][plotdf$Type == facet_labels[i]]
                }
                #only show for outliers
                plotdf$upper_lim[plotdf$upper_lim - plotdf$susie_effect < 5] = 0
                plotdf$lower_lim[plotdf$susie_effect - plotdf$lower_lim < 5] = 0
        }
        if(is.null(sd_names)){
        pdf(paste0("sim_rst/",filename))
        print(
        plotdf %>%
                ggplot(aes(x = true_effect, y = susie_effect)) +
                geom_point() +
                geom_abline(intercept = 0, slope = 1, color = "red") +
                #coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
                facet_wrap(~Type,scales = "free",nrow=1) +
                xlab("True gene-to-trait effect size") +
                ylab("SuSiE regression coefficient estimate") +
                theme_bw() +
                theme(text = element_text(size = 10, face = "bold"),aspect.ratio=1)
        )
        dev.off()
        }
        if(!is.null(sd_names)){
        pdf(paste0("sim_rst/",filename))
        print(
        plotdf %>%
                ggplot(aes(x = true_effect, y = susie_effect)) +
                geom_point() +
                geom_errorbar(aes(ymin=lower_lim, ymax=upper_lim)) +
                geom_abline(intercept = 0, slope = 1, color = "red") +
                #coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
                facet_wrap(~Type,scales = "free",nrow=1) +
                xlab("True gene-to-trait effect size") +
                ylab("SuSiE regression coefficient estimate") +
                theme_bw() +
                theme(text = element_text(size = 10, face = "bold"),aspect.ratio=1)
        )
        dev.off()
        }
}


plot_effects(c("susie_pip_pgene_grcp_dap_25_dap",
               "susie_pip_pgene_grcp_dap_25_elastic_net",
               "susie_pip_pgene_grcp_dap_25_smr"),
             c("susie_effect_pgene_grcp_dap_25_dap",
               "susie_effect_pgene_grcp_dap_25_elastic_net",
               "susie_effect_pgene_grcp_dap_25_smr"),
             c("ISuSiE",
               "PrediXcan",
               "SMR"),
             "isusie_effect_plots_vary_pred_model.pdf",
             c("susie_effect_sd_pgene_grcp_dap_25_dap",
               "susie_effect_sd_pgene_grcp_dap_25_elastic_net",
               "susie_effect_sd_pgene_grcp_dap_25_smr"))

plot_effects3(c("susie_pip_pgene_grcp_pi1_ptwas",
                "susie_pip_pgene_grcp_dap_25_dap",
              #  "susie_pip_pgene_grcp_dap_50_dap",
                "susie_pip_pgene_grcp_dap_75_dap",
                "susie_pip_pgene_grcp_dap_95_dap"),
              c("susie_effect_pgene_grcp_pi1_ptwas",
                "susie_effect_pgene_grcp_dap_25_dap",
              #  "susie_effect_pgene_grcp_dap_50_dap",
                "susie_effect_pgene_grcp_dap_75_dap",
                "susie_effect_pgene_grcp_dap_95_dap"),
             c("CPIP Threshold = 0",
               "CPIP Threshold = 0.25",
             #  "CPIP Threshold = 0.50",
               "CPIP Threshold = 0.75",
               "CPIP Threshold = 0.95"),
	      "isusie_effect_plots_0_95.pdf",
             c("susie_effect_sd_pgene_grcp_pi1_ptwas",
                "susie_effect_sd_pgene_grcp_dap_25_dap",
              #  "susie_effect_sd_pgene_grcp_dap_50_dap",
                "susie_effect_sd_pgene_grcp_dap_75_dap",
                "susie_effect_sd_pgene_grcp_dap_95_dap"))

cov <- c("25","50","75","95")

#SuSiE fine-mapping

plot_effects2(paste0("susie_pip_pgene_grcp_",cov,"_susie"),
              paste0("susie_effect_pgene_grcp_",cov,"_susie"),
              c("SuSiE coverage = 0.25",
               "SuSiE coverage = 0.50",
               "SuSiE coverage = 0.75",
               "SuSiE coverage = 0.95"),
              "isusie_susie_effects.pdf",
              paste0("susie_effect_sd_pgene_grcp_",cov,"_susie"))


#Uncensored DAP

plot_effects2(c("susie_pip_pgene_grcp_uncensored_25_dap",
                "susie_pip_pgene_grcp_uncensored_50_dap",
                "susie_pip_pgene_grcp_uncensored_75_dap",
                "susie_pip_pgene_grcp_uncensored_95_dap"),
              c("susie_effect_pgene_grcp_uncensored_25_dap",
                "susie_effect_pgene_grcp_uncensored_50_dap",
                "susie_effect_pgene_grcp_uncensored_75_dap",
               "susie_effect_pgene_grcp_uncensored_95_dap"),
             c("DAP CPIP threshold = 0.25",
               "DAP CPIP threshold = 0.50",
               "DAP CPIP threshold = 0.75",
               "DAP CPIP threshold = 0.95"),
             "dap_uncensored_effects.pdf",
             c("susie_effect_sd_pgene_grcp_uncensored_25_dap",
               "susie_effect_sd_pgene_grcp_uncensored_50_dap",
               "susie_effect_sd_pgene_grcp_uncensored_75_dap",
               "susie_effect_sd_pgene_grcp_uncensored_95_dap"))

#Uncensored SuSiE


plot_effects2(c("susie_pip_pgene_grcp_uncensored_25_susie",
                "susie_pip_pgene_grcp_uncensored_50_susie",
                "susie_pip_pgene_grcp_uncensored_75_susie",
                "susie_pip_pgene_grcp_uncensored_95_susie"),
              c("susie_effect_pgene_grcp_uncensored_25_susie",
                "susie_effect_pgene_grcp_uncensored_50_susie",
                "susie_effect_pgene_grcp_uncensored_75_susie",
               "susie_effect_pgene_grcp_uncensored_95_susie"),
             c("SuSiE coverage = 0.25",
               "SuSiE coverage = 0.50",
               "SuSiE coverage = 0.75",
               "SuSiE coverage = 0.95"),
             "susie_uncensored_effects.pdf",
             c("susie_effect_sd_pgene_grcp_uncensored_25_susie",
               "susie_effect_sd_pgene_grcp_uncensored_50_susie",
               "susie_effect_sd_pgene_grcp_uncensored_75_susie",
               "susie_effect_sd_pgene_grcp_uncensored_95_susie"))


#Individual plots (uncensored, PrediXcan, SMR, ISuSiE)

pdf("sim_rst/uncensored_individ.pdf")
d %>%
	mutate(upper = susie_effect_pgene_grcp_pi1_ptwas + 1.96 * susie_effect_sd_pgene_grcp_pi1_ptwas) %>%
	mutate(lower = susie_effect_pgene_grcp_pi1_ptwas - 1.96 * susie_effect_sd_pgene_grcp_pi1_ptwas) %>%
	mutate(upper = case_when(upper - susie_effect_pgene_grcp_pi1_ptwas < 5 ~ 0,TRUE ~ upper)) %>%
	mutate(lower = case_when(susie_effect_pgene_grcp_pi1_ptwas - lower < 5 ~ 0,TRUE ~ lower)) %>%
	ggplot(aes(x = true_effect,y = susie_effect_pgene_grcp_pi1_ptwas)) + 
	geom_point() +
 	geom_errorbar(aes(ymin=lower, ymax=upper)) +
	geom_abline(intercept = 0, slope = 1, color = "red") +
	xlab("True gene-to-trait effect size") +
	ylab("SuSiE regression coefficient estimate") +
	ggtitle("PTWAS") +
	theme_bw() +
	theme(text = element_text(size = 10, face = "bold"),aspect.ratio=1)
dev.off()

pdf("sim_rst/predixcan_individ.pdf")
d %>%
        ggplot(aes(x = true_effect,y = susie_effect_pgene_grcp_dap_25_elastic_net)) +
        geom_point() +
        geom_abline(intercept = 0, slope = 1, color = "red") +
        xlab("True gene-to-trait effect size") +
        ylab("SuSiE regression coefficient estimate") +
        ggtitle("PrediXcan") +
        theme_bw() +
        theme(text = element_text(size = 10, face = "bold"),aspect.ratio=1)
dev.off()

pdf("sim_rst/smr_individ.pdf")
d %>%
        ggplot(aes(x = true_effect,y = susie_effect_pgene_grcp_dap_25_smr)) +
        geom_point() +
        geom_abline(intercept = 0, slope = 1, color = "red") +
        xlab("True gene-to-trait effect size") +
        ylab("SuSiE regression coefficient estimate") +
        ggtitle("SMR") +
        theme_bw() +
        theme(text = element_text(size = 10, face = "bold"),aspect.ratio=1)
dev.off()

pdf("sim_rst/isusie_individ.pdf")
d %>%
        ggplot(aes(x = true_effect,y = susie_effect_pgene_grcp_dap_25_dap)) +
        geom_point() +
        geom_abline(intercept = 0, slope = 1, color = "red") +
        xlab("True gene-to-trait effect size") +
        ylab("SuSiE regression coefficient estimate") +
        ggtitle("ISuSiE") +
        theme_bw() +
        theme(text = element_text(size = 10, face = "bold"),aspect.ratio=1)
dev.off()


