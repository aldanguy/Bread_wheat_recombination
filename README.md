# scripts

Pipeline to estimate and compare LD-based recombination landscapes of 4 bread wheat landraces populations


How to run this pipeline ?

Create the following folder architecture
1) Create "scripts" folder and paste scripts of this Github repository


2) Create "amont" folder, which should contain:
- 161010_Chinese_Spring_v1.0_pseudomolecules_AGP.tsv (scaffold boundaries on RefSeq v1.0, reference: helene.rimbert@inrae.fr)
- csre_genotyping.txt (CsRe genotyping matrix, reference: pierre.sourdille@inrae.fr)
- BW_261K_wp3.ped and BW_261K_wp3.map (genotyping data of landraces, reference: sophie.bouchet@inrae.fr)
- accessions_4403_alice.csv (passeport data on landraces, reference: sophie.bouchet@inrae.fr)
- landrace632_8741hap.dis (pairwise simple matching distance matrix of landraces, reference: sophie.bouchet@inrae.fr)
- landrace632_8741hap_sophie.var (haplotypes matrix of landraces, reference: sophie.bouchet@inrae.fr)
- stat_haplo8741.txt (frontiers of haplotypic blocks, reference: sophie.bouchet@inrae.fr)
- Codes_chr.txt (ID of chromosomes in letters or numbers, reference: sophie.bouchet@inrae.fr)
- Decoupage_chr_ble.tab (frontiers of genomic regions, reference: sophie.bouchet@inrae.fr)
- SNP_positions.txt (positions of SNP on RefSeq v1.0, reference: sophie.bouchet@inrae.fr)
- domestication.csv (frontiers of domestication and agronomic genes, derived from Pont & al. 2019)
- iwgsc_refseqv1.0_HighConf_UTR_2017May05.gff3 (frontiers of genes, 5' and 3' UTR ... from URGI wbesite)

You can find all these raw files on ....



3) Create analysis folder (example with the French calendar date : "020820"). Should be empty.

The pipeline will keep using the bash script "base_020820.sh" to find appropriate names of folder and files during analysis. In this bash script, you should specify three variables: 
- r_amont should be equal to the path of your "amont" folder. Example : r_amont=/work/adanguy/these/pipeline/amont/
- r_scripts should be equal to the path of your "scripts" folder. Example r_scripts=/work/adanguy/these/pipeline/scripts/
- r should be equal to the path of your analysis folder. Example r=/work/adanguy/these/pipeline/020820/

Once you modified "base_020820.sh", you just have to run the bash script "first_020820.sh". You can use following code:
- source /work/adanguy/these/pipeline/scripts/base_020820.sh
- sbatch ${first}

This pipeline was designed to work on the remote servor "genologin.toulouse.inra.fr" (http://bioinfo.genotoul.fr/) using SLURM. 
It will call jobs arrays (up to 300 running at the same time) and need varying level memory (no more than 50 Gb for a specific job). 
If you want to run the complete pipeline, it will recquire ~ 40K hours of CPU running (around 34K hours only for PHASE estimation including 17K for specific SNP dataset and a little bit less for common SNP dataset).
Outputs are heavy (few hundreds of Gb).

Required installed softwares are:
- R (R-3.6.2 in all scripts except modele.sh using R-3.4.3_bis because of asreml R library)
- ASREML (R library)
- PLINK
- PHASE (path is /save/servin/bin/PHASE in PHASE_2.sh)
- HAPFLK (path is /save/servin/Envs/hapflkdev/bin/activate in FST.sh)
- bedtools-2.27.1 and bedops-v2.4.35 to use bedr library in HR_permutation.sh
