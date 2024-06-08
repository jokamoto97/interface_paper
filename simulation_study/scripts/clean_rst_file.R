library(data.table)
library(dplyr)

d <- fread("sim_data/sim.summary.pred.models.susie.p12.cis.pgwas.js.pgene.95.coloc.grcp.compare.coverage17.out",h=T)

keep <- c("sim_num","gene","true_effect","ptwas_z","elastic_net_z","smr_z","GLCP","GRCP","susie_pip_pgene_grcp_dap_25_dap","susie_pip_multi_gene_twas_dap_25","susie_effect_multi_gene_twas_dap_25","susie_pip_default_dap_25","susie_effect_default_dap_25","susie_pip_focus_dap_25","susie_effect_focus_dap_25","susie_pip_ctwas_priors_ptwas","susie_effect_ctwas_priors_ptwas","susie_pip_pgene_grcp_dap_25_dap","susie_effect_pgene_grcp_dap_25_dap","susie_pip_true_priors_censored_wts_ptwas","susie_effect_true_priors_censored_wts_ptwas","susie_pip_pgene_grcp_dap_25_dap","susie_effect_pgene_grcp_dap_25_dap","susie_pip_pgene_grcp_dap_25_elastic_net","susie_effect_pgene_grcp_dap_25_elastic_net","susie_pip_pgene_grcp_dap_25_smr","susie_effect_pgene_grcp_dap_25_smr","susie_pip_pgene_grcp_dap_25_dap","susie_effect_pgene_grcp_dap_25_dap","susie_pip_pgene_grcp_dap_50_dap","susie_effect_pgene_grcp_dap_50_dap","susie_pip_pgene_grcp_dap_75_dap","susie_effect_pgene_grcp_dap_75_dap","susie_pip_pgene_grcp_dap_95_dap","susie_effect_pgene_grcp_dap_95_dap","susie_pip_pgene_grcp_25_susie","susie_effect_pgene_grcp_25_susie","susie_pip_pgene_grcp_50_susie","susie_effect_pgene_grcp_50_susie","susie_pip_pgene_grcp_75_susie","susie_effect_pgene_grcp_75_susie","susie_pip_pgene_grcp_95_susie","susie_effect_pgene_grcp_95_susie","susie_pip_pgene_grcp_pi1_ptwas","susie_effect_pgene_grcp_pi1_ptwas","susie_pip_pgene_grcp_dap_25_dap","susie_effect_pgene_grcp_dap_25_dap","susie_pip_pgene_grcp_dap_50_dap","susie_effect_pgene_grcp_dap_50_dap","susie_pip_pgene_grcp_dap_75_dap","susie_effect_pgene_grcp_dap_75_dap","susie_pip_pgene_grcp_dap_95_dap","susie_effect_pgene_grcp_dap_95_dap",
          "susie_pip_pgene_grcp_dap_25_dap",
               "susie_pip_pgene_grcp_dap_25_elastic_net",
               "susie_pip_pgene_grcp_dap_25_smr","susie_effect_pgene_grcp_dap_25_dap",
               "susie_effect_pgene_grcp_dap_25_elastic_net",
               "susie_effect_pgene_grcp_dap_25_smr","susie_effect_sd_pgene_grcp_dap_25_dap",
               "susie_effect_sd_pgene_grcp_dap_25_elastic_net",
               "susie_effect_sd_pgene_grcp_dap_25_smr","susie_pip_pgene_grcp_pi1_ptwas",
                "susie_pip_pgene_grcp_dap_25_dap",
                "susie_pip_pgene_grcp_dap_50_dap",
                "susie_pip_pgene_grcp_dap_75_dap",
                "susie_pip_pgene_grcp_dap_95_dap","susie_effect_pgene_grcp_pi1_ptwas",
                "susie_effect_pgene_grcp_dap_25_dap",
                "susie_effect_pgene_grcp_dap_50_dap",
                "susie_effect_pgene_grcp_dap_75_dap",
                "susie_effect_pgene_grcp_dap_95_dap","susie_effect_sd_pgene_grcp_pi1_ptwas",
                "susie_effect_sd_pgene_grcp_dap_25_dap",
                "susie_effect_sd_pgene_grcp_dap_50_dap",
                "susie_effect_sd_pgene_grcp_dap_75_dap",
                "susie_effect_sd_pgene_grcp_dap_95_dap","susie_pip_pgene_grcp_uncensored_25_dap",
                "susie_pip_pgene_grcp_uncensored_50_dap",
                "susie_pip_pgene_grcp_uncensored_75_dap",
                "susie_pip_pgene_grcp_uncensored_95_dap","susie_effect_pgene_grcp_uncensored_25_dap",
                "susie_effect_pgene_grcp_uncensored_50_dap",
                "susie_effect_pgene_grcp_uncensored_75_dap",
               "susie_effect_pgene_grcp_uncensored_95_dap","susie_effect_sd_pgene_grcp_uncensored_25_dap",
               "susie_effect_sd_pgene_grcp_uncensored_50_dap",
               "susie_effect_sd_pgene_grcp_uncensored_75_dap",
               "susie_effect_sd_pgene_grcp_uncensored_95_dap","susie_pip_pgene_grcp_uncensored_25_susie",
                "susie_pip_pgene_grcp_uncensored_50_susie",
                "susie_pip_pgene_grcp_uncensored_75_susie",
                "susie_pip_pgene_grcp_uncensored_95_susie","susie_effect_pgene_grcp_uncensored_25_susie",
                "susie_effect_pgene_grcp_uncensored_50_susie",
                "susie_effect_pgene_grcp_uncensored_75_susie",
               "susie_effect_pgene_grcp_uncensored_95_susie","susie_effect_sd_pgene_grcp_uncensored_25_susie",
               "susie_effect_sd_pgene_grcp_uncensored_50_susie",
               "susie_effect_sd_pgene_grcp_uncensored_75_susie",
               "susie_effect_sd_pgene_grcp_uncensored_95_susie")


cov <- c("25","50","75","95")

keep <- c(keep,paste0("susie_pip_pgene_grcp_",cov,"_susie"),
              paste0("susie_effect_pgene_grcp_",cov,"_susie"),paste0("susie_effect_sd_pgene_grcp_",cov,"_susie"))

keep <- unique(keep)

out <- d %>% dplyr::select(all_of(keep))

write.table(out,"sim_data/rst_summary.txt",row.names=F,col.names=T,sep='\t',quote=F)  
