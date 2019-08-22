*************************************************************************
***           Scripts to perform the human analyses                   ***
***   from the paper Transition from background selection to          ***
***  associative overdominance promotes diversity in regions of       ***
***                low recombination                                  ***
***   by Gilbert, Pouyet, Excoffier and Peischl                       ***
***   Theses scripts were prepared by Fanny Pouyet on 2019.           ***
***              Last update 22.08.2019                               ***
*************************************************************************


#-------------------------- PI COMPUTATION  -----------------------------#
*This part relies on the following files: 
	1 .  	The vcf from 1000G project, 1 per chromosome (http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/).
	2 . 	The bed file that combine strictMack infos (ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/),
		the YRI recombination map from HapMap (ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/) and
		the CpG sites identified in the paper Pouyet et al., 2018 (referenced in the current manuscript; any C->T followed by a G and G->A following a C 
			from the genotype table of diallelic SNPs of Pouyet et al,2018).
		I used bedtools merge for strictMask and the rec. map and bedtools subtract for the merged bed and the CpG sites.
		I extracted using awk a bed file for regions with recombination rate > 0cM/Mb and <=0.05 cM/Mb as well as a bed for regions >= 1 and < 1.5 cM/Mb.
		--> Theses files are: strict_mask_YRI.recomb.Rec0p05-0p1.noCpG.bed and strict_mask_YRI.recomb.Rec1p-1p5.noCpG.bed
	3 .	For each population, I have the name of the 10 individuals studied in Pouyet et al., 2018 (files named: 1000G_${pop}names.txt
	4.	I used bedtools intersect to get the 13,385,820 SNPs from Pouyet et al., 2018 that are within strictMask and transorm them in 1 based coordinates : 1000GPAN.position.strict_mask_1based 
	
# I used the following commandline in a bash script to get pi per sites:

#!/bin/bash

file="nucleotide_diversity_strict_mask_sites_chr"
sites="1000GPAN.position.strict_mask_1based"
#i is the chromosome
i=$1
for pop in YRI LWK GBR IBS KHV JPT BEB CLM PJL PEL; do   
vcftools --gzvcf /storage/data/1000G/vcf/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --site-pi --positions ${sites} --keep 1000G_${pop}names.txt --out ${file}${pop}${i}
done


# Merge all these data to have a final table with PI for each population and all autosomes.

for chrom in {1..22};do 
 for pop in YRI LWK GBR IBS BEB CLM KHV JPT PJL PEL; do 
  cut -f3 nucleotide_diversity_strict_mask_sites_chr${pop}${chrom}.sites.pi > chr${pop}${chrom}.sites.pi
  done
 awk '{print $1"\t"$2"\t"$2+1}' nucleotide_diversity_strict_mask_sites_chrYRI${chrom}.sites.pi > chr${chrom}.pos
 paste chr${chrom}.pos chrYRI${chrom}.sites.pi chrLWK${chrom}.sites.pi chrGBR${chrom}.sites.pi chrIBS${chrom}.sites.pi chrBEB${chrom}.sites.pi chrPJL${chrom}.sites.pi chrKHV${chrom}.sites.pi chrJPT${chrom}.sites.pi chrCLM${chrom}.sites.pi chrPEL${chrom}.sites.pi >  chr${chrom}.sites.pi
 rm chr${chrom}.pos chrYRI${chrom}.sites.pi chrLWK${chrom}.sites.pi chrGBR${chrom}.sites.pi chrIBS${chrom}.sites.pi chrKHV${chrom}.sites.pi chrJPT${chrom}.sites.pi chrBEB${chrom}.sites.pi chrCLM${chrom}.sites.pi chrPJL${chrom}.sites.pi chrPEL${chrom}.sites.pi 
 sed -i '1d' chr${chrom}.sites.pi
 bedtools intersect -a chr${chrom}.sites.pi -b strict_mask_YRI.recomb.lowRec0p05.noCpG.bed> chr${chrom}.lowRec0p05.sites.pi 
 bedtools intersect -a chr${chrom}.sites.pi -b strict_mask_YRI.recomb.Rec1p-1p5.noCpG.bed > chr${chrom}.RR1p-1p5.sites.pi
 cat header chr${chrom}.lowRec0p05.sites.pi > tmp ; mv tmp chr${chrom}.lowRec0p05.sites.pi
 cat header chr${chrom}.RR1p-1p5.sites.pi > tmp ; mv tmp chr${chrom}.RR1p-1p5.sites.pi
done

*Then run windowingpi.R to have windows of 2Mb (step 0.5Mb)
#FINAL FILES :
	1. nucdiv.Rec0p05-0p1.windows2Mb-step0p5Mb.pi
	2. nucdiv.Rec0p05-1p0.windows2Mb-step0p5Mb.pi


#----------------------- MAKING Figures  --------------------#
*See the R script Fig_mAOD_HumanAnalyses.R
	1 .	Makes Supplemental Index with figures of Pi scans, SFS and heatmaps
	 	while Fig4 is a panel of supplemental item done using inkscape
	 	I have computed the number of SNPs for LR and MR for each window on the file NbSNP_LR_MR_perwindow.txt using the genotype table with all SNPs

