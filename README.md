# Scripts and computer code: Evolution of recombination landscapes in diverging populations of bread wheat

Pipeline to estimate and compare LD-based recombination landscapes of 4 bread wheat landraces populations


How to run this pipeline ?

Create the following folder architecture

1) Create "scripts" folder and paste scripts of this Github repository


2) Create "amont" folder, which should contain:

- CsRe genotyping matrix: csre_genotyping.txt (download at https://doi.org/10.5281/zenodo.4486612)

- Genotyping data of landraces and SNP positions on RefSeq v1.0 (download https://doi.org/10.5281/zenodo.4518374). : 
  . Collection_genotyping_matrix_410k.map
  . Collection_genotyping_matrix_410k.ped
  . SNP_positions_410k.txt 
  Data are under embargo until February 2022. Instead, you can use a subset of genotyping matrix available in Balfourier & al. 2019 supplementary material
  (https://advances.sciencemag.org/content/5/5/eaav0536, Data file S2. "Genotyping data of 4506 wheat accession with 113,457 genome-wide SNPs")

- Several supplementary files from Balfourier & al. 2019 (reference: sophie.bouchet @inrae.fr)
  . passeport data of landraces: accessions_4403_alice.csv
  . pairwise simple matching distance matrix of landraces: landrace632_8741hap.dis
  . haplotypes matrix of landraces: landrace632_8741hap_sophie.var
  . frontiers of haplotypic blocks: stat_haplo8741.txt

- Files about bread wheat genome (reference : helene.rimbert@inrae.fr)
  . ID of chromosomes in letters or numbers: Codes_chr.txt
  . frontiers of genomic regions: Decoupage_chr_ble.tab
  . frontiers of genes features: iwgsc_refseqv1.0_HighConf_UTR_2017May05.gff3
  . scaffold boundaries on RefSeq v1.0: 161010_Chinese_Spring_v1.0_pseudomolecules_AGP.tsv 
  
 - Info derived from Pont et al. 2019 about position of domestication and improvment genes: domestication.csv (https://doi.org/10.1038/s41588-019-0393-z)

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
- R (R-3.6.2 in all scripts except modele.sh using R-3.4.3_bis because of ASReml R library)
- ASReml (R library)
- PLINK
- PHASE (path is /save/servin/bin/PHASE in PHASE_2.sh)
- HAPFLK (path is /save/servin/Envs/hapflkdev/bin/activate in FST.sh)
- bedtools-2.27.1 and bedops-v2.4.35 to use bedr library in HR_permutation.sh

4) Main outputs are available in Zendo at https://doi.org/10.5281/zenodo.4486586


