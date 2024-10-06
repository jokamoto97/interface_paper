#Compare PCGs implicated by KBA, INTACT, and INTERFACE

library(data.table)
library(dplyr)
library(eulerr)
library(stringr)
library(tidyr)

#Limit to gene-metabolite pairs tested in cgene fine-mapping

region_files <- list.files("data/region_cid_genes")

get_cid_gene <- function(x){

        cid <- str_split_fixed(x,pattern = "_",n = 2)[1,1]

        region <- str_sub(str_split_fixed(x,pattern = "_",n = 2)[1,2],end = -5)

        #if(region == "6_22475112_79683053"){
        #       return(NULL)
        #}else{

        rst <- fread(paste0("data/region_cid_genes/",x))

        genes <- rst$GENES

        genes <- str_split_fixed(genes,pattern = "_",n = 2)[,1]

        return(paste0(cid,"_",genes))
        #}
}

test_pairs <- do.call(c,lapply(region_files,get_cid_gene))

#INTACT rst

intact_rst <- fread("data/intact_rst_sig05_t5.txt")

intact_pairs <- unique(paste0(str_split_fixed(intact_rst$Gene,pattern = "_",n = 3)[,1],"_",
                       str_split_fixed(intact_rst$Gene,pattern = "_",n = 3)[,2]))

#Cgene fine-mapping rst

cgene_fm <- fread("data/interface_rst_sig.txt")

fm_pairs <- paste0(cgene_fm$CID,"_",
                       str_split_fixed(cgene_fm$Gene,pattern = "_",n = 2)[,1])

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
        mutate(cid_gene = paste0(CID,"_",gene)) %>%
        filter(cid_gene %in% test_pairs)

kba_pairs <- unique(eric_dat$cid_gene)



#Make Venn Diagram

intact_fm <- intersect(intact_pairs, fm_pairs)

intact_kba <- intersect(intact_pairs, kba_pairs)

fm_kba <- intersect(fm_pairs, kba_pairs)

##Diagram entries

intact_only <- length(setdiff(intact_pairs, union(fm_pairs,kba_pairs)))

fm_only <- length(setdiff(fm_pairs, union(intact_pairs,kba_pairs)))

kba_only <- length(setdiff(kba_pairs, union(fm_pairs,intact_pairs)))

intact_fm_kba <- length(intersect(intact_fm, kba_pairs))

intact_and_fm <- length(setdiff(intact_fm,intersect(intact_fm,kba_pairs)))

intact_and_kba <- length(setdiff(intact_kba,intersect(intact_fm,kba_pairs)))

fm_and_kba <- length(setdiff(fm_kba,intersect(intact_fm,kba_pairs)))

length(union(kba_pairs,union(intact_pairs,fm_pairs)))

sum(c(intact_only,fm_only,kba_only,intact_and_fm,intact_and_kba,fm_and_kba,intact_fm_kba))


#All three approaches

pdf('results/compare_intact_fm_kba.pdf')
s2 <- c("INTACT" = intact_only,
        "INTERFACE" = fm_only,
        "Knowledge-Based" = kba_only,
        "INTACT&INTERFACE" = intact_and_fm,
        "INTACT&Knowledge-Based" = intact_and_kba,
        "INTERFACE&Knowledge-Based" = fm_and_kba,
        "INTACT&INTERFACE&Knowledge-Based" = intact_fm_kba)

plot(euler(s2),quantities = list(type = "counts",cex = 2),fills = list(fill = c( "red","steelblue4","grey"), alpha = 0.5),
     labels = list(col = "black",cex = 1.25))
dev.off()



#INTERFACE-only pairs

#metab_dat <- fread("data/annotationAll_FINAL_w_HMDB_ID.csv") %>% 
#	mutate(CID = paste0("C",CHEM_ID)) %>%
#	dplyr::select(CID,BIOCHEMICAL_NAME,TYPE,SUPER_PATHWAY,SUB_PATHWAY)

#data.frame("pair" = setdiff(fm_pairs, union(intact_pairs,kba_pairs))) %>%
#	separate(pair, into = c("CID", "Gene"),sep = '_') %>%
#	merge(metab_dat,by="CID",all.x=T) %>%
#	dplyr::select(Gene,CID,BIOCHEMICAL_NAME,TYPE,SUPER_PATHWAY,SUB_PATHWAY) %>%
#	#filter(!(startsWith(CID,"C9"))) %>%
#	write.table(file="interface_kba_difference_set.txt",sep='\t',row.names=F,col.names=T,quote=F)
