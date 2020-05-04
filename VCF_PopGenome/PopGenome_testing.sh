# testing how well I can analyze data using PopGenome

cd /cerberus/projects/chrwhe/PopGenome_testing

# files needed

# genome
ln -s /cerberus/projects/shared_napi_rapae/assemblies/Pieris_napi_fullAsm_chomoOnly.fasta .
# gtf
ln -s /cerberus/projects/shared_napi_rapae/Pieris_napi_annotation/gff/gene-build/rc3_pienap.sorted.gff .
# bed
cp /cerberus/projects/chrwhe/Pieris_napi_old_demography/bams/rg_files/chromosomes_2_25.bed .
cp /cerberus/projects/chrwhe/Pieris_napi_old_demography/bams/all.sorted_gt3.merged.sorted.readdepth_lt_4.N_repeat.merged.bed.gz .
# VCF file
cp /cerberus/projects/chrwhe/Pieris_napi_old_demography/bams/rg_files/all.freebayes.merged.Chromosome_2_25.SNPs.filtered.final.recode.vcf.gz .


# working with
https://cran.r-project.org/web/packages/PopGenome/vignettes/An_introduction_to_the_PopGenome_package.pdf


# https://github.com/tonig-evo/workshop-popgenome

cd /cerberus/projects/chrwhe/PopGenome_testing/testing_tutorial

git clone https://github.com/tonig-evo/workshop-popgenome.git

cd workshop-popgenome
unzip PopGenome_data.zip
