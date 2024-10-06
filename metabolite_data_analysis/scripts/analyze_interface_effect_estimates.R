#Analyze INTERFACE posterior mean effect sizes

library(data.table)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)

interface_rst <- fread("data/interface_rst.txt")

interface_rst$Gene_symbol <- str_split_fixed(interface_rst$Gene,pattern = "_",n=2)[,1]

interface_rst$cid_gene <- paste0(interface_rst$CID,"_",interface_rst$Gene_symbol)

interface_rst$KBA_pair <- FALSE

interface_rst$abs_effect <- abs(interface_rst$Posterior_Effect)

#Eric's KBA

eric_file <- 'data/tab_kba_genes_4_dap_signal_clusters.txt'
eric_dat <- fread(eric_file) %>%
        #filter(CID %in% intact_rst$metab_id) %>%
        dplyr::select(CHROM,START,STOP,CID,eric_gene_group)

tmp <- rbindlist(
  lapply(strsplit(eric_dat$eric_gene_group, "\\W"), function(x) data.table(t(x))),
  fill = TRUE
)

eric_dat <- cbind(eric_dat,tmp) %>%
        pivot_longer(cols = paste0("V",seq(1,9)),names_to = 'tmp',values_to = 'gene') %>%
        filter(is.na(gene) == F & gene != 'unknown') %>%
        mutate(cid_gene = paste0(CID,"_",gene))

kba_pairs <- unique(eric_dat$cid_gene)

interface_rst$KBA_pair[interface_rst$cid_gene %in% kba_pairs] <- TRUE





pdf("results/effect_size_distribution_kba.pdf")
interface_rst %>%
        mutate(KBA_pair = case_when(KBA_pair == TRUE ~ "KBA pair",
                                    KBA_pair == FALSE ~ "Not KBA pair")) %>%
        mutate(KBA_pair = factor(KBA_pair,levels = rev(c("KBA pair","Not KBA pair")))) %>%
        ggplot(aes(x = abs_effect,fill=KBA_pair)) +
        geom_histogram(binwidth=0.1) +
        xlab("Absolute INTERFACE gene-to-metabolite effect") +
        ylab("Number of gene-metabolite pairs") +
        coord_cartesian(xlim = c(0,4),ylim = c(0,200)) +
        theme_bw() +
        theme(text = element_text(size = 10,face="bold"),aspect.ratio=1,legend.title=element_blank())
dev.off()


t.test(interface_rst$abs_effect[interface_rst$KBA_pair == T],interface_rst$abs_effect[interface_rst$KBA_pair == F])

