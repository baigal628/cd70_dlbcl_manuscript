#!/bin/bash
#SBATCH --job-name=encypt_data
#SBATCH --mem=60G        # total memory need
#SBATCH -n 32 #Number of cores
#SBATCH -N 1 # on one node
#SBATCH --error=/liulab/galib/hodgkin_lymphoma_mtang/EGA_data/lsf_%j_%x.err
#SBATCH --output=/liulab/galib/hodgkin_lymphoma_mtang/EGA_data/lsf_%j_%x.err
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=galib@ds.dfci.harvard.edu # email to notify



java -jar /liulab/galib/hodgkin_lymphoma_mtang/hodgkin_lymphoma_cellranger_output/EGA-Cryptor-2.0.0/ega-cryptor-2.0.0.jar -i /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/all_analysis.tar.gz -o /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/encryped/ -t 64

