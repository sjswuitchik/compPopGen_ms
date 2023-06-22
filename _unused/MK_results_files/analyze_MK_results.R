##preliminary analysis of SNPIRE/MK data for PEQG poster

library(tidyverse)

birds<-tibble(spcode=c("anaPla", "cotJap", "egrGar", "falPer", "nipNip", "taeGut"), species_name=c("Anas platyrhynchos", "Coturnix japonica", "Egretta garzetta", "Falco peregrinus", "Nipponia nippon", "Taeniopygia guttata"), common_name=c("Mallard", "Japanese quail", "Little egret", "Peregrine falcon", "Crested ibis", "Zebra finch"))
fishes<-tibble(spcode=c("angAng", "cluHar", "cypCar", "funHet", "hapBur", "hipCom", "larCro", "punPun"), species_name=c("Anguilla anguilla", "Clupea harengus", "Cyprinus carpio", "Fundulus heteroclitus", "Haplochromis burtoni", "Hippocampus comes", "Larimichthys crocea", "Pungitius pungitius"), common_name=c("European eel", "Atlantic herring", "Eurasian carp", "Mummichog", "Burton's mouthbrooder", "Tiger tail seahorse", "Large yellow croaker", "Ninespine stickleback"))
reptiles<-tibble(spcode=c("anoCar"), species_name=c("Anolis carolinensis"), common_name=c("Green anole"))

all_species <- bind_rows(birds, fishes, reptiles)

#read in all files

read_snipre<-function(spcode, path) {
  file<-paste0(path, "/", spcode, "_snipre_output.tsv")
  read_tsv(file, col_names=TRUE) %>% mutate(spcode = spcode) %>% mutate(gene = toupper(gene))
}

watterson<-function(K, n, bp) {
  an <- sum(1 / 1:(n-1))
  return(K/an/bp)
}

path<-"/Users/tim/Projects/popgen/compPopGen_ms/MK_results_files"

bird_results <- lapply(birds$spcode, read_snipre, path) %>% bind_rows()
fish_results <- lapply(fishes$spcode, read_snipre, path) %>% bind_rows()
rep_results <- lapply(reptiles$spcode, read_snipre, path) %>% bind_rows()
all_res <- bind_rows(bird_results, fish_results, rep_results)


#some basic overviews

all_res %>% group_by(spcode) %>% summarize(syn_div = sum(FS), tot_syn=sum(Tsil)) %>% full_join(all_species) %>% mutate(br = syn_div / tot_syn)

sum_stats <- all_res %>% 
  filter(total_poly > 0) %>% 
  group_by(spcode) %>% 
  summarize(syn_div = sum(FS),
            syn_poly = sum(PS),
            tot_syn=sum(Tsil),
            num_pos_snipre = sum(SnIPRE.class == "pos"),
            num_pos_mk = sum(mk_pval < 0.05 & dos > 0),
            tot_genes = length(gene),
            nsamp = round(median(npop))) %>% 
  mutate(br = syn_div / tot_syn,
         prop_pos_mk = num_pos_mk / tot_genes,
         prop_pos_s = num_pos_snipre / tot_genes,
         an = sapply(nsamp, function(x) sum(1/1:(x-1))),
         theta = syn_poly / an / tot_syn) %>%
  full_join(all_species)


filter(all_res, total_poly > 0, total_div > 0) %>% ggplot(aes(x=SnIPRE.est, y=dos)) + geom_point()

filter(all_res, spcode == "cotJap", total_poly > 0) %>% ggplot(aes(x=SnIPRE.est, y=alpha, color=SnIPRE.class)) + geom_smooth()
filter(all_res, spcode == "angAng", total_poly > 0) %>% ggplot(aes(x=SnIPRE.est, y=alpha, color=SnIPRE.class)) + geom_smooth()
filter(all_res, spcode == "taeGut", total_poly > 0) %>% ggplot(aes(x=SnIPRE.est, y=alpha)) + geom_smooth()

filter(all_res, spcode %in% c("angAng", "taeGut", "cotJap", "egrGar", "hapBur", "larCro"), total_poly > 0, total_div > 0) %>% ggplot(aes(x=SnIPRE.est, y=dos, color=spcode)) + geom_smooth()


#compare to old results
mk_res_comp %>% filter(species %in% c("Egarzetta", "Nnippon", "Tguttata"), PS+PR>0) %>% mutate(dos = FR/(FR + FS) - PR/(PR + PS)) %>% ggplot(aes(y=dos, x=SnIPRE.est, color=species)) + geom_smooth()

all_res %>% filter(spcode %in% c("egrGar", "nipNip", "taeGut")) %>% mutate(dosnew = FR/(FR + FS) - PR/(PR + PS)) %>% ggplot(aes(y=dosnew, x=SnIPRE.est, color=spcode)) + geom_smooth()

all_res %>% filter(spcode == "egrGar", total_poly > 0, total_div > 0) %>% ggplot(aes(x=dos)) + geom_histogram(bins=100)
mk_res_comp %>% filter(species == "Egarzetta", PR+PS > 0, FR+FS > 0) %>% mutate(dos = FR/(FR + FS) - PR/(PR + PS)) %>% ggplot(aes(x=dos)) + geom_histogram(bins=100)

all_res %>% filter(spcode == "nipNip", total_poly > 0, total_div > 0) %>% ggplot(aes(x=dos)) + geom_histogram(bins=100)
mk_res_comp %>% filter(species == "Nnippon", PR+PS > 0, FR+FS > 0) %>% mutate(dos = FR/(FR + FS) - PR/(PR + PS)) %>% ggplot(aes(x=dos)) + geom_histogram(bins=100)

taeGut_comp <- full_join(filter(all_res, spcode == "taeGut"), filter(mk_res_comp,species=="Tguttata") %>% mutate(dos = FR/(FR + FS) - PR/(PR + PS)), by=c("gene" = "gene_name"), suffix=c(".new", ".old"))
taeGut_comp %>% ggplot(aes(x=dos.old, y=dos.new)) + geom_density_2d()
taeGut_comp %>% ggplot(aes(x=SnIPRE.est.old, y=SnIPRE.est.new)) + geom_density_2d()
taeGut_comp %>% mutate(div.old = FR.old+FS.old, div.new = FR.new+FS.new) %>% filter(div.old <3000, div.new<3000) %>% ggplot(aes(x=div.old, y=div.new)) + geom_point() + geom_abline(slope = 1, intercept = 0)
taeGut_comp %>% ggplot(aes(x=PR.old+PS.old, y=PR.new+PS.new)) + geom_point()  + geom_abline(slope = 1, intercept = 0)

egrGar_comp <- full_join(filter(all_res, spcode == "egrGar"), filter(mk_res_comp,species=="Egarzetta") %>% mutate(dos = FR/(FR + FS) - PR/(PR + PS)), by=c("gene" = "gene_name"), suffix=c(".new", ".old"))
egrGar_comp %>% ggplot(aes(x=dos.old, y=dos.new)) + geom_density_2d()
egrGar_comp %>% ggplot(aes(x=SnIPRE.est.old, y=SnIPRE.est.new)) + geom_density_2d()
egrGar_comp %>% mutate(div.old = FR.old+FS.old, div.new = FR.new+FS.new) %>% filter(div.old <3000, div.new<3000) %>% ggplot(aes(x=div.old, y=div.new)) + geom_point() + geom_abline(slope = 1, intercept = 0)
egrGar_comp %>% ggplot(aes(x=PR.old+PS.old, y=PR.new+PS.new)) + geom_point()  + geom_abline(slope = 1, intercept = 0)

nipNip_comp <- full_join(filter(all_res, spcode == "nipNip"), filter(mk_res_comp,species=="Nnippon") %>% mutate(dos = FR/(FR + FS) - PR/(PR + PS)), by=c("gene" = "gene_name"), suffix=c(".new", ".old"))
nipNip_comp %>% ggplot(aes(x=dos.old, y=dos.new)) + geom_density_2d()
nipNip_comp %>% ggplot(aes(x=SnIPRE.est.old, y=SnIPRE.est.new)) + geom_density_2d()
nipNip_comp %>% mutate(div.old = FR.old+FS.old, div.new = FR.new+FS.new) %>% filter(div.old <3000, div.new<3000) %>% ggplot(aes(x=div.old, y=div.new)) + geom_point() + geom_abline(slope = 1, intercept = 0)
nipNip_comp %>% ggplot(aes(x=PR.old+PS.old, y=PR.new+PS.new)) + geom_point()  + geom_abline(slope = 1, intercept = 0)

cotJap_comp <- full_join(filter(all_res, spcode == "cotJap"), filter(mk_res_comp,species=="Cjaponica") %>% mutate(dos = FR/(FR + FS) - PR/(PR + PS)), by=c("gene" = "gene_name"), suffix=c(".new", ".old"))
cotJap_comp %>% ggplot(aes(x=dos.old, y=dos.new)) + geom_density_2d()
cotJap_comp %>% ggplot(aes(x=SnIPRE.est.old, y=SnIPRE.est.new)) + geom_density_2d()
cotJap_comp %>% ggplot(aes(x=PR.old+PS.old, y=PR.new+PS.new)) + geom_point()  + geom_abline(slope = 1, intercept = 0)
cotJap_comp %>% mutate(div.old = FR.old+FS.old, div.new = FR.new+FS.new) %>% filter(div.old <3000, div.new<3000) %>% ggplot(aes(x=div.old, y=div.new)) + geom_point() + geom_abline(slope = 1, intercept = 0)


#figure out complete genes

bird_genes <- bird_results %>% count(gene) %>% rename(bird_n = n)
fish_genes <- fish_results %>% count(gene) %>% rename(fish_n = n)
rep_genes <- rep_results %>% count(gene) %>% rename(reptile_n = n)

genes <- full_join(bird_genes, fish_genes) %>% full_join(rep_genes) 
genelist <- filter(genes, bird_n == 6) %>% pull(gene)

#initial analysis set
all_res_filt_class.1 <- all_res %>% filter(gene %in% genelist, spcode %in% birds$spcode) %>% 
  group_by(gene) %>% count(SnIPRE.class) %>% pivot_wider(names_from = SnIPRE.class, values_from = n, values_fill = 0) %>%
  mutate(prop_pos = pos / 6)


all_res_filt_class.1 <- all_res %>% filter(gene %in% genelist, spcode %in% birds$spcode) %>% 
  group_by(gene) %>% mutate(pos_num = sum(mk_pval<0.05)) %>% select(gene, pos_num) %>% distinct()


#do a little bit of GO testing on old vs new

##taeGut##

taeGut.new <- all_res %>% filter(spcode == "taeGut") %>% select(gene, FR, FS, PR, PS, dos, alpha, total_div, total_poly, mk_pval, mk_pval_fdr)
taeGut.old <- mk_res_comp %>% filter(species == "Tguttata") %>% 
  mutate(dos = FR/(FR + FS) - PR/(PR + PS),
         total_div = FR + FS,
         total_poly = PR + PS) %>% 
  select(gene_name, FR, FS, PR, PS, dos, alpha, total_div, total_poly, mk_pval, mk_pval_fdr)

#go enrichment

library(DOSE)
library(clusterProfiler)
library(biomaRt)
library(org.Gg.eg.db)

fore_id.new <- taeGut.new %>% filter(mk_pval < 0.05) %>% pull(gene) %>% bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb = org.Gg.eg.db)
uni_id.new <- taeGut.new %>% pull(gene) %>% bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb = org.Gg.eg.db)


fore_id.old <- taeGut.old %>% filter(mk_pval < 0.05) %>% pull(gene_name) %>% bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb = org.Gg.eg.db)
uni_id.old <- taeGut.old %>% pull(gene_name) %>% bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb = org.Gg.eg.db)

taeGut_kegg.new <- enrichKEGG(fore_id.new$ENTREZID,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=.2,universe=uni_id.new$ENTREZID,keyType="ncbi-geneid")

taeGut_kegg.old <- enrichKEGG(fore_id.old$ENTREZID,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=.2,universe=uni_id.old$ENTREZID,keyType="ncbi-geneid")

enrichGO(gene=fore_id.old$ENTREZID, OrgDb=org.Gg.eg.db, universe = uni_id.old$ENTREZID)
enrichGO(gene=fore_id.new$ENTREZID, OrgDb=org.Gg.eg.db, universe = uni_id.$ENTREZID)

##double check numbers 

mk_old.poster <- mk_res_comp %>% 
  filter(species != "Gvarius") %>%
  group_by(gene_name) %>% count(SnIPRE.class) %>% 
  pivot_wider(names_from = "SnIPRE.class", values_from = "n", values_fill = 0) %>%
  mutate(n = neut+pos+neg) %>% filter(n == 5)


#calculate expectations from random permutations

pa<-mk_res_comp %>% filter(species == "Cjaponica") %>% summarize(pos = sum(SnIPRE.class == "pos")/n()) %>% pull(pos)
pb<-mk_res_comp %>% filter(species == "Egarzetta") %>% summarize(pos = sum(SnIPRE.class == "pos")/n()) %>% pull(pos)
pc<-mk_res_comp %>% filter(species == "Falbicollis") %>% summarize(pos = sum(SnIPRE.class == "pos")/n()) %>% pull(pos)
pd<-mk_res_comp %>% filter(species == "Nnippon") %>% summarize(pos = sum(SnIPRE.class == "pos")/n()) %>% pull(pos)
pe<-mk_res_comp %>% filter(species == "Tguttata") %>% summarize(pos = sum(SnIPRE.class == "pos")/n()) %>% pull(pos)

#make simulated data
mk_perms<-tibble(perm = rep(1:1000, each=7596), 
                 spa=as.vector(replicate(1000,rbinom(7596,1,pa))),
                 spb=as.vector(replicate(1000,rbinom(7596,1,pb))),
                 spc=as.vector(replicate(1000,rbinom(7596,1,pc))),
                 spd=as.vector(replicate(1000,rbinom(7596,1,pd))),
                 spe=as.vector(replicate(1000,rbinom(7596,1,pe))))


mk_perms.res <- mk_perms %>% mutate(pos = spa+spb+spc+spd+spe) %>% group_by(perm) %>% count(pos) %>% pivot_wider(names_from = "pos", values_from = "n", values_fill = 0, names_prefix = "pos")

table(mk_old.poster$pos)



conf_limits <- tibble(num = c(1,2,3,4),
                      median = c(median(mk_perms.res$pos1),
                                 median(mk_perms.res$pos2),
                                 median(mk_perms.res$pos3),
                                 median(mk_perms.res$pos4)),
                      lowerb = c(quantile(mk_perms.res$pos1, probs = 0.025),
                                 quantile(mk_perms.res$pos2, probs = 0.025),
                                 quantile(mk_perms.res$pos3, probs = 0.025),
                                 quantile(mk_perms.res$pos4, probs = 0.025)),
                      upperb = c(quantile(mk_perms.res$pos1, probs = 0.975),
                                 quantile(mk_perms.res$pos2, probs = 0.975),
                                 quantile(mk_perms.res$pos3, probs = 0.975),
                                 quantile(mk_perms.res$pos4, probs = 0.975))) %>% filter(num > 1)
                      
real_data <- tibble(num = c(1,2,3,4), obs = as.numeric(table(mk_old.poster$pos)[2:5])) %>% filter(num > 1)                   

ggplot(real_data, aes(x=num, y=obs)) + geom_point(size=6, color="purple") +
  geom_line(data=conf_limits, aes(x=num, y=median), size=1, linetype="dotted") +
  geom_line(data=conf_limits, aes(x=num, y=lowerb), linetype="dashed", color="blue", size=1.8) +  
  geom_line(data=conf_limits, aes(x=num, y=upperb), linetype="dashed", color="blue", size=1.8) +
  theme_classic(base_size = 30) + xlab("Number of species") +
  ylab("Number of genes") + scale_x_continuous(breaks=c(2,3,4))

fore<-mk_old.poster %>% filter(pos >= 2) %>% pull(gene_name)
back<-mk_old.poster %>% pull(gene_name)

library(DOSE)
library(clusterProfiler)
library(org.Gg.eg.db)

fore.ncbi<-bitr(fore,fromType="SYMBOL", toType="ENTREZID", OrgDb = org.Gg.eg.db)
back.ncbi<-bitr(back,fromType="SYMBOL", toType="ENTREZID", OrgDb = org.Gg.eg.db)

enrich<-enrichKEGG(fore.ncbi$ENTREZID,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=.2,universe=back.ncbi$ENTREZID,keyType="ncbi-geneid")
dotplot(enrich)
