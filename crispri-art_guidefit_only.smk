import glob
import os
import pandas as pd

# stationary parameters for now
PROJECT="/groups/doudna/projects/BA_Cas_Bioinformatics/CRISPRi-ARTv1"
RUN="20240516_CRISPRiART-418868809__MS2only"
ADAPTERS="/home/badler/bin/CRISPRi-ART/params/default_adapter.fasta"
FAST2Q__US="GGTTTGAAAC"
FAST2Q__DS="ATGCTTGGGC"
LT = 0.1
UT = 0.9

metadata_file = "{project}/data/metadata/{run}__metadata.tsv".format(project=PROJECT,run=RUN)
df_metadata = pd.read_table(metadata_file)
sample_list = df_metadata["sample_name"].tolist()
lib_list = list(df_metadata["lib_name"].unique())
extended_sample_list = [sample[1].lib_name + '/' + sample[1].sample_name for sample in df_metadata.iterrows()]

#inelegant solution for file renaming before snakemake. Shouldn't matter if already done.
for dir_name in os.listdir("{project}/data/reads/{run}".format(project=PROJECT,run=RUN)) : 
    if any(sample_name in dir_name for sample_name in sample_list):
        new_dir_name = dir_name.split("_L1")[0] + "_L1"
        os.rename("{project}/data/reads/{run}/{dir_name}".format(project=PROJECT,run=RUN,dir_name=dir_name),"{project}/data/reads/{run}/{new_dir_name}".format(project=PROJECT,run=RUN,new_dir_name=new_dir_name))

rule all:
    input: 
        #count reads
        expand("{project}/data/tmp/{run}/counting_completed.txt".format(project=PROJECT,run=RUN)),
        #process, normalize, and determine guide fitness from reads
        expand("{project}/data/fitness/{run}/{lib_name}/fitness_checkpoint.txt", project=PROJECT,run=RUN,lib_name=lib_list),

rule trim_reads_wrapper: 
    input:
        # trim reads
        run = expand("{project}/data/merged_reads/{run}/{extended_sample_name}_merged.fastq", project=PROJECT,run=RUN,extended_sample_name=extended_sample_list)
    output: 
        completed = "{project}/data/tmp/{run}/trimming_completed.txt".format(project=PROJECT,run=RUN)
    params:
        base_dir = directory(expand("{project}/data/tmp".format(project=PROJECT))),
        base_dir_run = directory(expand("{project}/data/tmp/{run}".format(project=PROJECT,run=RUN)))
    shell:
        """
        mkdir -p {params.base_dir}
        mkdir -p {params.base_dir_run}
        touch {output.completed}
        """

rule trim_reads:
    input:
        r1 = lambda wildcards: glob.glob("{project}/data/reads/{run}/{sample_name}_L1/{sample_name}*_L001_R1_001.fastq.gz".format(project=PROJECT,run=RUN,sample_name=wildcards.extended_sample_name.split('/')[-1])),
        r2 = lambda wildcards: glob.glob("{project}/data/reads/{run}/{sample_name}_L1/{sample_name}*_L001_R2_001.fastq.gz".format(project=PROJECT,run=RUN,sample_name=wildcards.extended_sample_name.split('/')[-1]))
    output: 
        merged = "{project}/data/merged_reads/{run}/{extended_sample_name}_merged.fastq",
        html_report = "{project}/data/merged_reads/{run}/reports/{extended_sample_name}_fastp.html",
        json_report = "{project}/data/merged_reads/{run}/reports/{extended_sample_name}_fastp.json",
    params:
        adapters = "{adapters}".format(adapters=ADAPTERS),
        base_dir = directory(expand("{project}/data/merged_reads".format(project=PROJECT))),
        base_dir_run = directory(expand("{project}/data/merged_reads/{run}".format(project=PROJECT,run=RUN)))
    shell: 
        """
        mkdir -p {params.base_dir}
        mkdir -p {params.base_dir_run}
        fastp \
            --in1 {input.r2} \
            --in2 {input.r1} \
            --merge \
            --merged_out {output.merged} \
            --adapter_fasta {params.adapters}\
            --trim_poly_g \
            --html {output.html_report} \
            --json {output.json_report}
        """

rule count_reads_wrapper: 
    input:
        checkpoint = "{project}/data/tmp/{run}/trimming_completed.txt".format(project=PROJECT,run=RUN),
        run = expand("{project}/data/counts/{run}/{lib_name}/fast2q_output/{lib_name}_fast2q_counts.csv", project=PROJECT,run=RUN,lib_name=lib_list)
    output: 
        completed = "{project}/data/tmp/{run}/counting_completed.txt".format(project=PROJECT,run=RUN)
    params:
        base_dir = directory(expand("{project}/data/tmp".format(project=PROJECT))),
        base_dir_run = directory(expand("{project}/data/tmp/{run}".format(project=PROJECT,run=RUN)))
    shell:
        """
        mkdir -p {params.base_dir}
        mkdir -p {params.base_dir_run}
        touch {output.completed}
        """

rule count_reads:
    output:
        fast2q_counts = "{project}/data/counts/{run}/{lib_name}/fast2q_output/{lib_name}_fast2q_counts.csv",
        fast2q_counts_stats = "{project}/data/counts/{run}/{lib_name}/fast2q_output/{lib_name}_fast2q_counts_stats.csv"
    params: 
        merged_base = lambda wildcards: directory("{project}/data/merged_reads/{run}/{lib_name}".format(project=PROJECT,run=RUN,lib_name=wildcards.lib_name)),
        guide_lib = lambda wildcards: "{project}/oligos/{lib_name}_gRNAs.csv".format(project=PROJECT,lib_name=wildcards.lib_name),
        base_dir = directory("{project}/data/counts/".format(project=PROJECT)),
        base_dir_run = directory(expand("{project}/data/counts/{run}".format(project=PROJECT,run=RUN))),
        base_dir_run_lib = lambda wildcards: directory(expand("{project}/data/counts/{run}/{lib_name}".format(project=PROJECT,run=RUN,lib_name=wildcards.lib_name))),
        US = FAST2Q__US,
        DS = FAST2Q__DS,
        FN = lambda wildcards: "{lib_name}_fast2q_counts".format(lib_name=wildcards.lib_name)
    shell:
        """
        mkdir -p {params.base_dir}
        mkdir -p {params.base_dir_run}
        mkdir -p {params.base_dir_run_lib}
        python -m fast2q -c \
            --s {params.merged_base} \
            --g {params.guide_lib} \
            --o {params.base_dir_run_lib} \
            --fn {params.FN} \
            --l 20 \
            --m 2 \
            --l 30 \
            --us {params.US} \
            --ds {params.DS} \
            --ph 20
        most_recent_folder=$(ls -td "{params.base_dir_run_lib}"/*/ | head -1)
        echo $most_recent_folder
        mv $most_recent_folder/* {params.base_dir_run_lib}/fast2q_output/
        ls {params.base_dir_run_lib}/fast2q_output/
        rm -R $most_recent_folder/
        """
rule guide_fitness: 
    input: 
        reads_checkpoint = expand("{project}/data/tmp/{run}/counting_completed.txt".format(project=PROJECT,run=RUN)),
        fast2q_counts = lambda wildcards: "{project}/data/counts/{run}/{lib_name}/fast2q_output/{lib_name}_fast2q_counts.csv".format(project=PROJECT,run=RUN,lib_name=wildcards.lib_name)
    output: 
        guide_fitness_checkpoint = "{project}/data/fitness/{run}/{lib_name}/fitness_checkpoint.txt"
    params: 
        base_dir = directory(expand("{project}/data/fitness/".format(project=PROJECT))),
        base_dir_run = directory(expand("{project}/data/fitness/{run}".format(project=PROJECT,run=RUN))),
        base_dir_run_lib = lambda wildcards: directory(expand("{project}/data/fitness/{run}/{lib_name}".format(project=PROJECT,run=RUN,lib_name=wildcards.lib_name))),
        base_dir_run_lib_guide = lambda wildcards: directory(expand("{project}/data/fitness/{run}/{lib_name}/guide_fitness".format(project=PROJECT,run=RUN,lib_name=wildcards.lib_name))),
        metadata = "{project}/data/metadata/{run}__metadata.tsv".format(project=PROJECT,run=RUN),
        lower_threshold = LT,
        upper_threshold = UT
    shell:
        """
        mkdir -p {params.base_dir}
        mkdir -p {params.base_dir_run}
        mkdir -p {params.base_dir_run_lib}
        mkdir -p {params.base_dir_run_lib_guide}
        python crispriart_guide_fitness__only.py \
            --counts {input.fast2q_counts} \
            --metadata {params.metadata} \
            --output {params.base_dir_run_lib} \
            --lower {params.lower_threshold} \
            --upper {params.upper_threshold} 
        """

        





