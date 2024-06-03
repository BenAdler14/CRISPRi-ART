# ---------------------------------------------------------------------------- #
# Packages and function import
# ---------------------------------------------------------------------------- #

import os
import shutil
import argparse
#from dna_features_viewer import GraphicFeature
#from dna_features_viewer import GraphicRecord
#from dna_features_viewer import BiopythonTranslator
import pandas as pd
#from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------- #
# Input arguments (CHECK ME - ESPECIALLY PARSER)
# ---------------------------------------------------------------------------- #
def get_args():
    parser = argparse.ArgumentParser(
        description=""" 
        # expected inputs
        # args.lib_name
        # args.project (optional, but would otherwise need to include all path input arguments manually)
        # args.run (optional, but would otherwise need to include all path input arguments manually)
        
        # path input arguments (optional, but would otherwise need to include project, run arguments)
        # args.input_metadata_file
        # args.fast2q_counts (optional, can generate from p
        """
    )

    # Required arguments
    parser.add_argument(
        "-i",
        "--input_fitness",
        type=str,
        required=True,
        default="",
        help="""
        Path to the input gene_fitness file.
        """,
    )

    parser.add_argument(
        "-c",
        "--condition",
        type=str,
        required=True,
        default="",
        help="""
        Which condition you are analyzing.
        """,
    )
    parser.add_argument(
        "-f",
        "--fitness_threshold",
        type=str,
        required=True,
        default="",
        help="""
        Which gene fitness threshold you want to use for determining Fit status.
        """,
    )

    parser.add_argument(
        "-p",
        "--pval",
        type=str,
        required=True,
        default="",
        help="""
        Which p_value threshold you want to use for KS Statistic.
        """,
    )
        
    # Outputs - only one is required.    
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        default="",
        help="""
        Name of the output directory.
        """,
    )

    args = parser.parse_args()

    # Validate inputs
        
    if (args.input_fitness == ""
    ):
        msg = "Please specify a input gene_fitness file_path!  "
        raise ValueError(msg)
        
    if (args.condition == ""
    ):
        msg = "Please specify a condition to analyze!  "
        raise ValueError(msg)
        
    if (args.fitness_threshold == ""
    ):
        msg = "Please specify a gene fitness threshold you want to use!  "
        raise ValueError(msg)

    if (args.pval == ""
    ):
        msg = "Please specify a p-value threshold you want to use!  "
        raise ValueError(msg)

    if (args.output == ""
    ):
        msg = "Please specify an output directory you want to use!  "
        raise ValueError(msg)

    return args

# ---------------------------------------------------------------------------- #
# Importing functions
# ---------------------------------------------------------------------------- #
def import_gene_fitness(path_gene_fitness) : 
    """
    Import gene fitness table. 
    """
    df_gene_fitness = pd.read_table(path_gene_fitness)
    return df_gene_fitness

def generate_analysis_id(condition, fitness_threshold, pval_threshold) : 
    """
    generate a file_id unique to this analysis
    """
    return "{condition}__fitness_{fitness_threshold}__pval_{pval_threshold}".format(condition=condition,fitness_threshold=fitness_threshold,pval_threshold=pval_threshold)

def initialize_output(path_output_dir, condition, fitness_threshold, pval_threshold) : 
    """
    Set up output directory to accommodate all output files. Return the base
    output path. 
    """
    if not os.path.exists(path_output_dir):
        os.makedirs(path_output_dir)
    
    # generate analysis-specific output directory - will overwrite an existing directory if it has been tried before
    analysis_id = generate_analysis_id(condition, fitness_threshold, pval_threshold)
    base_output_path = os.path.join(path_output_dir,analysis_id)

    if os.path.exists(base_output_path):
        shutil.rmtree(base_output_path)
    os.makedirs(base_output_path)
    return base_output_path
    
    
    
# ---------------------------------------------------------------------------- #
# Fitness thresholding functions
# ---------------------------------------------------------------------------- #
def filter_gene_fitness(df_gene_fitness, condition, fitness_threshold, pval_threshold, base_output_path) :
    """
    Filter df_gene_fitness by the user-provided thresholds (upper only) and output
    a modified gene fitness file including "Fit", "Semi-Fit" and "Not Fit" values. 
    """
    default_columns = ['gene', 'gRNA_name', 'gRNA_start', 'gRNA_end', 'num_guides', condition, condition+'_KSpvalue']
    df_gene_fitness_FIT = df_gene_fitness[default_columns].copy(deep=True)
    df_gene_fitness_FIT['Fit genes'] = 'Not Fit'
    
    # select fit_genes and save table as tsv for fitness - use for essentiality labels elsewhere = UPDATED
    df_gene_fitness_FIT.loc[(df_gene_fitness_FIT[condition] > float(fitness_threshold)) & (df_gene_fitness_FIT[condition+'_KSpvalue'] <= float(pval_threshold)),'Fit genes'] = 'Fit'
    df_gene_fitness_FIT.loc[(df_gene_fitness_FIT[condition] > float(fitness_threshold)) & (df_gene_fitness_FIT[condition+'_KSpvalue'] > float(pval_threshold)),'Fit genes'] = 'Semi-Fit'

    analysis_id = generate_analysis_id(condition, fitness_threshold, pval_threshold)
    path_output_table = os.path.join(base_output_path, analysis_id+'.tsv')
    df_gene_fitness_FIT.to_csv(path_output_table,index=False,sep='\t')
    
    return df_gene_fitness_FIT

# ---------------------------------------------------------------------------- #
# Plotting functions
# ---------------------------------------------------------------------------- #
def volcano_plot_v4(df_gene_fitness_FIT, condition, fitness_threshold, pval_threshold, base_output_path) :
    analysis_id = generate_analysis_id(condition, fitness_threshold, pval_threshold)
    path_plot_out = os.path.join(base_output_path, analysis_id + '_volcano_plot_v4.svg')

    # -log10 transform pvalues
    df_gene_fitness_FIT[condition+'_KSpvalue'] = -1*np.log10(df_gene_fitness_FIT[condition+'_KSpvalue'])
    
    # Select data for highlighting gene fitness
    x_notfit = df_gene_fitness_FIT[df_gene_fitness_FIT['Fit genes'] == 'Not Fit'][condition].values
    y_notfit = df_gene_fitness_FIT[df_gene_fitness_FIT['Fit genes'] == 'Not Fit'][condition+'_KSpvalue'].values
    x_semifit = df_gene_fitness_FIT[df_gene_fitness_FIT['Fit genes'] == 'Semi-Fit'][condition].values
    y_semifit = df_gene_fitness_FIT[df_gene_fitness_FIT['Fit genes'] == 'Semi-Fit'][condition+'_KSpvalue'].values
    x_fit = df_gene_fitness_FIT[df_gene_fitness_FIT['Fit genes'] == 'Fit'][condition].values
    y_fit = df_gene_fitness_FIT[df_gene_fitness_FIT['Fit genes'] == 'Fit'][condition+'_KSpvalue'].values
    
    # Set up plotting
    plt.figure(figsize=(5,4))
    ax = plt.subplot()

    #first all the genes
    plt.plot(x_notfit, y_notfit, color = '#6EA6CD', alpha=0.5, marker = '.', linestyle = ' ')

    # color all the fit, but not significant genes
    plt.plot(x_semifit, y_semifit, color = '#F67E4B', alpha=0.5, marker = '.', linestyle = ' ')
    
    #color the all significant genes
    plt.plot(x_fit, y_fit, color = '#A50026', alpha=0.5, marker = '.', linestyle = ' ')

    # title and axis labels
    title = analysis_id
    ax.set_title(title, fontsize = 16)
    #ax.set_xlim(-4.25,4.25)
    #ax.set_ylim(-.5,12.5)
    plt.axhline((np.log10(float(pval_threshold))*-1), linestyle='--', color='black', linewidth = 1)
    plt.axvline(float(fitness_threshold), linestyle='--', color='black', linewidth = 1)
    ax.set_xlabel("$Log_{2}$ (Fold Change)", fontsize=16)
    ax.set_ylabel("$-Log_{10}$ (p-value)", fontsize=16)
    ax.spines['top'].set_visible(False) 
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.tick_params(axis = 'both', labelsize=14, bottom=True)
    plt.savefig(path_plot_out, format = 'svg')
    
    return

    
# ---------------------------------------------------------------------------- #
# Main (CHECK ME)
# ---------------------------------------------------------------------------- #
def main():
    # set fonts for plotting
    #plt.rcParams['font.sans-serif'] = "Arial"
    #plt.rcParams['font.family'] = "sans-serif"
    
    print("Reading inputs...")
    args = get_args()
    base_output_path = initialize_output(args.output, args.condition, args.fitness_threshold, args.pval)
    df_gene_fitness = import_gene_fitness(args.input_fitness)

    print("Determining fit genes...")
    # Separating genes by fitness
    df_gene_fitness_FIT = filter_gene_fitness(df_gene_fitness, args.condition, args.fitness_threshold, args.pval, base_output_path)

    print("Plotting fitness...")
    volcano_plot_v4(df_gene_fitness_FIT, args.condition, args.fitness_threshold, args.pval, base_output_path)

    print('Fit genes with fitness threshold of {fitness_threshold} and p-value of {pval}: '.format(fitness_threshold=args.fitness_threshold,pval=args.pval))
    print(len(df_gene_fitness_FIT[df_gene_fitness_FIT['Fit genes'] == 'Fit']))
    print()
    print('Semi-Fit genes with fitness threshold of {fitness_threshold} and p-value of {pval}: '.format(fitness_threshold=args.fitness_threshold,pval=args.pval))
    print(len(df_gene_fitness_FIT[df_gene_fitness_FIT['Fit genes'] == 'Semi-Fit']))
    print()
    print('Not Fit genes with fitness threshold of {fitness_threshold} and p-value of {pval}: '.format(fitness_threshold=args.fitness_threshold,pval=args.pval))
    print(len(df_gene_fitness_FIT[df_gene_fitness_FIT['Fit genes'] == 'Not Fit']))
    print()
    
if __name__ == "__main__":
    main()    