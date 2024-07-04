# CRISPRi-ART

This is a simple pipeline used to go from paired end reads to gene fitness for transcriptome-wide CRISPRi-ART screening.  

If you find CRISPRi-ART useful in your work, please cite: 

**Genome-wide Characterization of Diverse Bacteriophages Enabled by RNA-Binding CRISPRi**

Benjamin A. Adler*,  Muntathar J. Al-Shimary*, Jaymin R. Patel, Emily Armbruster,  David Colognori,  Emeric J. Charles, Kate V. Miller,  Arushi Lahiri,  Marena Trinidad,  Ron Boger,  Jason Nomburg,  Sebastien Beurnier,  Michael L. Cui, Rodolphe Barrangou,  Vivek K. Mutalik,  Joseph S. Schoeniger,  Joseph A. Pogliano,  David F. Savage,  Jennifer A. Doudna+,  Brady F. Cress+

# Installation 
1. Clone repo into your directory 
> git clone https://github.com/BenAdler14/CRISPRi-ART.git

2. create a conda environment to import CRISPRi-ART
> conda env create -f envs/CRISPRi-ART.yaml

# Usage 
1. Activate the CRISPRi-ART conda environment
> conda activate CRISPRi-ART

2. Create a project directory named "MY_PROJECT" with a sequencing run named "RUN". An example is shown for a project named "SAMPLE_PROJECT" for a sequencing run named "20240325_CRISPRiART_Requeue-725488425". The associated sample_sheet used for Illumina sequencing is shown under SAMPLE_PROJECT/data/sample_sheets/2024_03_25_Multi_lib_samplesheet.csv

3. Upload your demultiplexed reads to MY_PROJECT/data/reads/RUN/

4. Upload your metadata file to MY_PROJECT/data/metadata/RUN__metadata.tsv . An example can be found under SAMPLE_PROJECT/data/metadata/20240325_CRISPRiART_Requeue-725488425__metadata.tsv

5. Update the PROJECT, RUN, and ADAPTERS parameters in crispri-art_defaults.smk file.

6. To generate fitness scores:
> snakemake -s crispri-art_defaults.smk -cores 8

7. Upon completion, guide fitness output and gene should be found at (where LIBRARY could be a unique CRISPRi-ART library):
> SAMPLE_PROJECT/data/fitness/RUN/LIBRARY/guide_fitness/
> SAMPLE_PROJECT/data/fitness/RUN/LIBRARY/gene_fitness/

8. Data can be analyzed here without further processing, but if you wish to make Fit, Semi-Fit, or Not Fit fitness calls you can run the following where CONDITION is an experimental condition for further analysis (see examples in metadata), FIT_THRESH is your desired fitness threshold and PVAL_THRESH is your desired pvalue:
> python crispri-art_fitness_thresholding.py --input_fitness SAMPLE_PROJECT/data/fitness/RUN/LIBRARY/gene_fitness/fitness_genes.tsv --condition CONDITION --fitness_threshold FIT_THRESH --pval PVAL_THRESH --output SAMPLE_PROJECT/data/fitness/RUN/LIBRARY/gene_fitness_analysis/
