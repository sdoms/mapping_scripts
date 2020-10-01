#chrX
sbatch --array=1-5 slurm_glm_ChrX_DNA_class.slurm
sbatch --array=1-4 slurm_glm_ChrX_DNA_phylum.slurm
sbatch --array=1-6%3 slurm_glm_ChrX_DNA_order.slurm
sbatch --array=1-12%4 slurm_glm_ChrX_DNA_family.slurm
sbatch --array=1-28%4 slurm_glm_ChrX_DNA_genus.slurm
sbatch --array=1-67%4 slurm_glm_ChrX_DNA_otu.slurm

sbatch --array=1-67%4 slurm_glm_ChrX_RNA_otu.slurm
sbatch --array=1-28%4 slurm_glm_ChrX_RNA_genus.slurm
sbatch --array=1-12%4 slurm_glm_ChrX_RNA_family.slurm
sbatch --array=1-6%3 slurm_glm_ChrX_RNA_order.slurm
sbatch --array=1-4 slurm_glm_ChrX_RNA_phylum.slurm
sbatch --array=1-5 slurm_glm_ChrX_RNA_class.slurm


# RNA
sbatch --array=1-4 rna_phylum.slurm
sbatch --array=1-5%3 rna_class.slurm
sbatch --array=1-12%4 rna_family.slurm
sbatch --array=1-6%3 rna_order.slurm
sbatch --array=1-28%4 rna_genus.slurm
sbatch --array=1-67%4 rna_otu.slurm

# DNA
sbatch --array=1-4 dna_phylum.slurm
sbatch --array=1-5%3 dna_class.slurm
sbatch --array=1-12%4 dna_family.slurm
sbatch --array=1-6%3 dna_order.slurm
sbatch --array=1-28%4 dna_genus.slurm
 sbatch --array=1-67%4 dna_otu.slurm
# histology
#sbatch --array=9,13,14,15,16,20,21,22,23,24,25,26,27,29,32,35,38,41,44 glm_hist.slurm
#sbatch --array=9,13,14,15,16,20,21,22,23,24,25,26,27,29,32,35,38,41,44 glm_chrX_hist.slurm
sbatch --array=8,13,14,20,21 glm_hist.slurm
sbatch --array=8,13,14,20,21 glm_chrX_hist.slurm
