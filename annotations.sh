#!/bin/bash

base=${1}


source ${base}


### Proximity of HR with genes



job_out=${r_log}annotations_genes.out
    
job_name=genes
	    
job_genes=$(sbatch -o ${job_out} -J ${job_name} --mem=50G --parsable ${r_scripts}annotations_genes.sh ${base}) 
    
    







### Proximity of HR with functions


job_out=${r_log}annotations_go.out
    
job_name=go
	    
# job_go=$(sbatch -o ${job_out} -J ${job_name} --mem=50G --parsable ${r_scripts}annotations_go.sh ${base}) 
    
    

job_out=${r_log}annotations_te.out
    
job_name=te
	    
# job_te=$(sbatch -o ${job_out} -J ${job_name} --mem=50G --parsable ${r_scripts}annotations_te.sh ${base}) 
    
    


