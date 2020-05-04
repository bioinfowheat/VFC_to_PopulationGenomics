rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures
## set your working directory, e.g. setwd("~/")
setwd("/cerberus/projects/chrwhe/PopGenome_testing/testing_tutorial/workshop-popgenome/")



# https://github.com/tonig-evo/workshop-popgenome

cd /cerberus/projects/chrwhe/PopGenome_testing/testing_tutorial

git clone https://github.com/tonig-evo/workshop-popgenome.git

cd workshop-popgenome
unzip PopGenome_data.zip

#
library(PopGenome)
# Reading data, read fasta from folder
GENOME.class <- readData("fasta")
get.sum.data(GENOME.class)
get.individuals(GENOME.class)

# Available statistics and examples
show.slots(GENOME.class)
# Run necessary module
GENOME.class <- F_ST.stats(GENOME.class)
GENOME.class <- neutrality.stats(GENOME.class)
GENOME.class@n.sites
GENOME.class@Pi
GENOME.class@Tajima.D

# question3
GENOME.class <- linkage.stats(GENOME.class)
GENOME.class@Wall.B

# Question 5, for per site estimate of Pi
GENOME.class@Pi / GENOME.class@n.valid.sites

# Available region data and statistics
GENOME.class@region.data
GENOME.class@region.stats
# Examples
GENOME.class@region.data@biallelic.sites[[1]][1:10]
GENOME.class@region.data@transitions[[1]][1:10]

# Question 6, number of gapped sites
GENOME.class@region.data@sites.with.gaps
length(GENOME.class@region.data@sites.with.gaps[[1]]) # 1454

# Question 7, singletons
GENOME.class@region.data@n.singletons
sum(GENOME.class@region.data@n.singletons[[1]]) # 67 - note the [[]] produced in the previous command


# Without defining populations
get.individuals(GENOME.class)
GENOME.class <- neutrality.stats(GENOME.class,detail=TRUE)
get.neutrality(GENOME.class)[[1]]
# Define populations with lists
GENOME.class <- set.populations(GENOME.class,list(
  c("CON","KAS-1","RUB-1","PER-1","RI-0","MR-0","TUL-0"),
  c("MH-0","YO-0","ITA-0","CVI-0","COL-2","LA-0","NC-1") ))
# Check whether grouping is set correctly
GENOME.class@region.data@populations
GENOME.class@region.data@populations2
GENOME.class@region.data@outgroup
# Recalculate statistics for populations
GENOME.class <-neutrality.stats(GENOME.class,detail=TRUE)
GENOME.class@Tajima.D
# Each population
get.neutrality(GENOME.class)[[1]]
get.neutrality(GENOME.class)[[2]]
# Set an outgroup
GENOME.class <-set.outgroup(GENOME.class,c("Alyrata"))
GENOME.class@region.data@outgroup
GENOME.class <- neutrality.stats(GENOME.class,detail=TRUE)
get.neutrality(GENOME.class)[[1]]
get.neutrality(GENOME.class)[[2]]

10.What do you have to pay attention to when applying the McDonald-Kreitman test? (see Whole_genome_analyses_using_VCF_files.pdf)
The outgroup needs to be defined as a population, not as an outgroup

###########
#   VCF
##########
rm(list=ls()) #clears all variables
# To read a VCF file using readVCF it needs to be compressed with bgzip and indexed with tabix.
# The tabix files need to be placed in the same folder as the vcf file.
# What parameters need to be defined
GENOME2.class <-readVCF("great_tit/vcf/LGE22.vcf.gz", 6000,"chrLGE22_Parus_Major_build_1.0.2",1,773534)
# arguments
# filename	the corresponding tabixed VCF-file
# numcols	  number of SNPs that should be read in as a chunk
# tid	      which chromosome ? (character)
# frompos	  start of the region
# topos	    end of the region

GENOME2.class@region.names
GENOME2.class <- neutrality.stats(GENOME2.class, FAST=TRUE)
get.sum.data(GENOME2.class)
GENOME2.class@region.data


get.neutrality(GENOME2.class)[[1]]

GENOME2.class <- F_ST.stats(GENOME2.class)
GENOME2.class@theta_Watterson # 1161.065
GENOME2.class@Pi # 1153.292
GENOME2.class@Tajima.D # -2.529237
# per bp Pi
GENOME2.class@Pi / GENOME2.class@n.sites
# per bp Pi, but the above might also include sites for which we had no data

GENOME2.class@region.data


#########
rm(list=ls()) #clears all variables
GENOME2.class <- readData("great_tit/vcf2",format="VCF", gffpath="great_tit/gff")
get.sum.data(GENOME2.class)
GENOME2.class@region.data
GENOME2.class <- set.synnonsyn(GENOME2.class, ref.chr="great_tit/fasta/LGE22.fasta")
GENOME2.class@region.data@synonymous
GENOME2.class@region.data@CodingSNPS
GENOME2.class.syn <- neutrality.stats(GENOME2.class,subsites="syn")
GENOME2.class.syn@Tajima.D
GENOME2.class.syn@theta_Watterson
#
GENOME2.class.nonsyn <- neutrality.stats(GENOME2.class,subsites="nonsyn")
GENOME2.class.nonsyn@Tajima.D
GENOME2.class.nonsyn@theta_Watterson

#######
# for all genes, TajD
GENOME2.class.genes <- neutrality.stats(GENOME2.class, subsites="gene")
GENOME2.class.genes@Tajima.D
GENOME2.class.genes@theta_Watterson
# per gene TajD
GENOME2.class.by_genes <- splitting.data(GENOME2.class, subsites="gene")
GENOME2.class.by_genes <- neutrality.stats(GENOME2.class.by_genes)
GENOME2.class.by_genes@Tajima.D
GENOME2.class.by_genes@theta_Watterson
# per gene NS TajD
GENOME2.class.by_genes_non <- neutrality.stats(GENOME2.class.by_genes,subsites="nonsyn")
GENOME2.class.by_genes_non@Tajima.D
# plot
plot(GENOME2.class.by_genes@Tajima.D , ylim=c(0,1), xlab="genes", ylab="Tajima's D",pch=3)
plot(GENOME2.class.by_genes@theta_Watterson,  xlab="genes", ylab="theta_Watterson",pch=3)

# but now I want these results with gene names attached
gene.positions <- get_gff_info(gff.file="great_tit/gff/LGE22.gff", chr="LGE22", feature="gene")
# check
is(gene.positions)
GENOME.class.split <- splitting.data(GENOME2.class, positions=gene.positions, type = 2)
# set type=2 as gene.positions contains genomic positions
gene.ids <- get_gff_info(gff.file="great_tit/gff/LGE22.gff", chr="LGE22", extract.gene.names=TRUE )
GENOME.class.split <- neutrality.stats(GENOME.class.split)
GENOME.class.split@Tajima.D

# wrangling this into a table.
library(tidyverse)
#
outfile <- data.frame(GENOME.class.split@Tajima.D,gene.ids)  %>%
  separate(gene.ids, c("geneID","Dbxref","Name","rest"), ";", extra = "merge") %>%
  rename(TajD=pop.1) %>% select(TajD,geneID,Name)
#
write_delim(outfile,"outfile.tsv",delim="\t")
