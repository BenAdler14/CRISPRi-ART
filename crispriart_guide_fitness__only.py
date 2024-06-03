# ---------------------------------------------------------------------------- #
# Packages and function import
# ---------------------------------------------------------------------------- #
import os
import argparse
import re
import pandas as pd
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------- #
# Input arguments
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
    
    '''parser.add_argument(
        "-l",
        "--lib_name",
        type=str,
        required=None,
        default="",
        help="""
        Name of the library used.
        """,
    )'''

    parser.add_argument(
        "-c",
        "--counts",
        type=str,
        required=True,
        default="",
        help="""
        Path to the input fastq counts file.
        """,
    )

    parser.add_argument(
        "-m",
        "--metadata",
        type=str,
        required=True,
        default="",
        help="""
        Path to the input experimental metadata file.
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

    # Optional arguments.    
    parser.add_argument(
        "-lt",
        "--lower",
        type=str,
        required=False,
        default="0.05",
        help="""
        Lower percentile of reads to drop from analysis. Default: 0.05
        """,
    )

    parser.add_argument(
        "-ut",
        "--upper",
        type=str,
        required=False,
        default="0.95",
        help="""
        Lower percentile of reads to drop from analysis. Default: 0.95
        """,
    )

    args = parser.parse_args()

    # Validate inputs
    '''if (args.lib_name == ""
    ):
        msg = "Please specify a lib_name file_path!  "
        raise ValueError(msg)''';
        
    if (args.counts == ""
    ):
        msg = "Please specify a counts file_path!  "
        raise ValueError(msg)
        
    if (args.metadata == ""
    ):
        msg = "Please specify a metadata file_path!  "
        raise ValueError(msg)
        
    if (args.output == ""
    ):
        msg = "Please specify a output file_path!  "
        raise ValueError(msg)

    if (args.lower == "0.05"
    ):
        msg = "Removing gRNAs from fitness represented below 0.05 fraction of read counts in base condition (default)"

    if (args.upper == "0.95"
    ):
        msg = "Removing gRNAs from fitness represented above 0.95 fraction of read counts in base condition (default)"

    return args

# ---------------------------------------------------------------------------- #
# Importing functions
# ---------------------------------------------------------------------------- #
def import_counts(path_counts) : 
    """
    Import counts table. 
    """
    df_counts = pd.read_csv(path_counts)
    df_counts = df_counts.rename(columns={"#Feature":"gRNA_name"})
    return df_counts

def import_metadata(path_metadata) : 
    """
    Import metadata table and convert to a cond_dict. 
    """
    df_metadata = pd.read_table(path_metadata)
    return df_metadata

def filter__metadata(df_metadata,df_samples) :
    """
    filter metadata based on the following principles
    1. Must be of type "base" or "sample" (not time0 - this is not used as of 1/2023)
    2. Must be contained in df_samples (note: with read processing pipelines, often has extra stuff appended)
    """
    df_metadata = df_metadata[df_metadata["sample_type"].isin(["base","sample"])]
    df_metadata = df_metadata[df_metadata["sample_name"].apply(lambda x : any(re.search(x+r'*',sample_col) for sample_col in df_samples.columns.tolist()))]

    return df_metadata

def fetch_sample_values(df_samples, list_metadata_values) : 
    """
    fetch equivalent column names from a df_samples table given a list of metadata values (which are a strict subset of df_sample
    column names)
    """
    list_metadata_regex = [re.compile(metadata_value+r".*") for metadata_value in list_metadata_values]
    list_sample_values = []
    for pattern in list_metadata_regex : 
        for sample in df_samples.columns : 
            if pattern.match(sample) : 
                list_sample_values.append(sample)
    return list_sample_values

def metadata__2__cond_dict(df_metadata,df_samples) : 
    """
    Programmatically convert a metadata dataframe to a condition dictionary. 
    The keys of the dictionary are formatted as 
    dict{lib_name__sample_set : {"base":{base_group_name : list_base_columns},
                                "sample":{sample1_group_name : list_sample1_columns,
                                                          ....                     ,
                                          sampleN_group_name : list_sampleN_columns}
                                }
    
    lib_name should be a superset of sample set (that is all members of "sample_set" should employ the same lib_name). However, a lib_name can employ several sample_sets because the employ different base.
    group names are of format: {sample_type}--{varX_name}_{varX_value}__{varY_name}_{varY_value}... # for as many variables as present in the metadata
    
    to the left of the "--" delimeter is organizational for automated comparisions. 
    to the right of the "--" delimeter is naming for automated comparisons.
    """
    
    cond_dict = {}
    var_names = sorted([col.split('__')[0] for col in df_metadata.columns if '__name' in col])
    
    for sample_set_name, sample_set_group in df_metadata.groupby(["lib_name","sample_set"]) : 
        sample_set_key = '{lib_name}__{sample_set}'.format(lib_name=sample_set_group['lib_name'].values[0],sample_set=sample_set_group['sample_set'].values[0])
        cond_dict[sample_set_key] = {"base":{},"sample":{}}
        
        for name, group in sample_set_group.groupby([col for col in sample_set_group.columns if col != 'sample_name']) : 
            group_unique = group[[col for col in group.columns if col != 'sample_name']].drop_duplicates()
            sample_type = group_unique['sample_type'].values[0]
            group_key = ""
            
            for var in var_names : 
                group_key += '{varX__name}_{varX__value}__'.format(varX__name=group_unique[var+'__name'].values[0],varX__value=group_unique[var+'__value'].values[0])
            group_key = group_key[:-2]

            metadata_values = group['sample_name'].tolist()
            group_values = fetch_sample_values(df_samples, metadata_values)
            cond_dict[sample_set_key][sample_type][group_key] = group_values
    return cond_dict    

# ---------------------------------------------------------------------------- #
# Fitness calculation functions
# ---------------------------------------------------------------------------- #
def generate_fitness_outputs(path_out,sample_set): 
    dict_output = {}
    dict_output["counts_QC"] = os.path.join(path_out,sample_set+"__counts_QC.csv")
    dict_output["counts_norm"] = os.path.join(path_out,sample_set+"__counts_norm.tsv")
    dict_output["guide_fitness"] = os.path.join(path_out,sample_set+"__fitness_guides.tsv")
    return dict_output

def fit_calc(df_counts, cond_dict, lower_percentile, upper_percentile, path_out) : 
    """
    Wrapper function to calculate log2FC guide fitness. 
    1. Generate output paths -- currently this is prepared for multiple sample_set combinations (ie unique base values) despite currently only being one right now. 
    2. Mask erroneous reads by the base_condition (if any step here bears optimization it's this one)
    3. Normalize reads (by arithmetic mean)
    4. Determine guide fitness
    """
    sample_set_list = []
    for sample_set in cond_dict.keys() :
        # generate output paths
        dict_output = generate_fitness_outputs(path_out,sample_set)

        # list_all_samples
        list_all_samples = []
        list_base = list(cond_dict[sample_set]["base"].values())[0]
        list_samples = []
        for sample_cond in cond_dict[sample_set]["sample"].keys() : 
            list_samples.extend(cond_dict[sample_set]["sample"][sample_cond])
        list_all_samples.extend(list_base)
        list_all_samples.extend(list_samples)

        # only modify df_counts based on the sample_set
        final_cols = ["gRNA_name"]
        final_cols.extend(list_all_samples)
        df_counts_sample_set = df_counts[final_cols].copy()

        # mask extreme values for reads
        df_counts_sample_set = counts_extrema_mask(df_counts_sample_set, dict_output["counts_QC"], list_base, lower_percentile=lower_percentile, upper_percentile=upper_percentile, drop_criteria='any')

        # normalize all read counts
        df_counts_sample_set_norm = counts_norm_within_sample(df_counts_sample_set, list_all_samples, path_out_norm=dict_output["counts_norm"], style='mean-1e6')
        
        # calculate guide fitness
        df_guide_fit = log2FC_naive(df_counts_sample_set_norm, cond_dict[sample_set], path_out_guide_fit=dict_output["guide_fitness"])

        # append to sample_set_list
        sample_set_list.append(sample_set)
    
    return sample_set_list

def counts_extrema_mask(df_counts, path_out_masked, list_base, lower_percentile=0.05, upper_percentile=1, drop_criteria='any') : 
    # for denominator columns, determine low/upper thresholds
    low_thresh = df_counts[list_base].quantile(lower_percentile)
    high_thresh = df_counts[list_base].quantile(upper_percentile)
    
    # create new df_counts that meets filtering criteria and merged by drop criteria
    keep_index = df_counts[list_base][(df_counts[list_base] > low_thresh) & (df_counts[list_base] < high_thresh)].dropna(how='any').index
    
    df_counts['QC_pass'] = [False]*len(df_counts)
    df_counts.loc[keep_index,'QC_pass'] = True
    
    df_counts.to_csv(path_out_masked,index=False)
    
    return df_counts


# normalize counts per index to 1e6 pseudocounts after adding 1 count to all guides. Other normalization options below. 
def counts_norm_within_sample(df_counts, list_cond_full, path_out_norm='', style='mean-1e6') : 
    """
    apply a normalization method. 
    default to mean-1e6 because it is the most straightforward normalization to read_depth method. 
    """
    
    df_counts[list_cond_full] += 1
    if style == 'mean-1e6' : 
        df_counts[list_cond_full] = df_counts[list_cond_full].apply(lambda x: 1e6 * x/x.sum(), axis=0)
    elif style == 'mean' : 
        df_counts[list_cond_full] = df_counts[list_cond_full].apply(lambda x: x/x.sum(), axis=0)
    elif style == 'minmax-1e6' : 
        df_counts[list_cond_full] = df_counts[list_cond_full].apply(lambda x: 1e6 * (x-x.min())/(x.max()-x.min()), axis=0)
    elif style == 'minmax' : 
        df_counts[list_cond_full] = df_counts[list_cond_full].apply(lambda x: (x-x.min())/(x.max()-x.min()), axis=0)
    elif style == 'geomean' :
        df_counts[list_cond_full] = df_counts[list_cond_full].apply(lambda x: x/gmean(x), axis=0)
    elif style == 'median-1e6' :
        df_counts[list_cond_full] = df_counts[list_cond_full].apply(lambda x: 1e6 * x/x.median(), axis=0)
    elif style == 'median' :
        df_counts[list_cond_full] = df_counts[list_cond_full].apply(lambda x: x/x.median(), axis=0)
    elif style == 'median ratio' : 
        df_counts[list_cond_full] = df_counts[list_cond_full].apply(lambda x: x/(np.median(x/gmean(x))), axis=0)        
    else:
        print('Not a valid within-sample normalization option')
    
    if path_out_norm != '' : 
        df_counts.to_csv(path_out_norm,sep='\t',index=False)
        
    return df_counts


def log2FC_naive(df_counts_norm, cond_dict_sample_set, path_out_guide_fit='', norm_con="base") :
    """
    Quick calculation of log2FC for all samples within a sample_set (unique base condition). Normalize to columns in the base condition.
    """
    
    # change below line depending on the final file structure
    df_guide_fit = df_counts_norm.copy()
    cond_base = list(cond_dict_sample_set["base"].keys())[0]
    col_base = cond_dict_sample_set["base"][cond_base]
    df_base = df_guide_fit[col_base]
    log2_base = np.mean(np.log2(df_base),axis=1)

    # get all samples
    list_cond_base = list(cond_dict_sample_set["base"].keys())
    list_cond_samples = []
    list_all_samples = []
    list_base = cond_dict_sample_set["base"][cond_base]
    list_all_samples.extend(list_base)
    for sample in list(cond_dict_sample_set["sample"].keys()) :
        list_cond_samples.append(sample)
        list_all_samples.extend(cond_dict_sample_set["sample"][sample])
    
    # calculate log2fc against the base
    df_guide_fit[list_all_samples] = df_guide_fit[list_all_samples].apply(lambda x: np.log2(x) - log2_base, axis=0)

    # calculate mean values for the base condition and t-stat like pval (should be near zero with near-one pval)
    df_base = df_guide_fit[col_base]
    df_cond = df_guide_fit[list_base]
    df_guide_fit[list_cond_base[0]] = np.mean(df_cond,axis=1)
    df_guide_fit[list_cond_base[0] + '__pval'] = stats.ttest_ind(df_base, df_cond, alternative='two-sided',axis=1).pvalue

    # for each non-base condition calculate mean values for the base condition and t-stat like pval
    for cond_group in list_cond_samples : 
        df_cond = df_guide_fit[cond_dict_sample_set["sample"][cond_group]]
        df_guide_fit[cond_group] = np.mean(df_cond,axis=1)
        df_guide_fit[cond_group + '__pval'] = stats.ttest_ind(df_base, df_cond, alternative='two-sided',axis=1).pvalue
    
    if path_out_guide_fit != '' : 
        df_guide_fit.to_csv(path_out_guide_fit,sep='\t',index=False)
        
    return df_guide_fit
    
# ---------------------------------------------------------------------------- #
# Write checkpoint
# ---------------------------------------------------------------------------- #
def write_checkpoint(sample_set_list, path_out) : 
    with open(os.path.join(path_out,"guide_fitness.txt"),"w") as f : 
        for sample_set in sample_set_list : 
            f.write(sample_set + '\n')
    f.close()
    return
    
# ---------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------- #
def main():
    args = get_args()

    print("Reading inputs...")
    df_counts = import_counts(args.counts)
    df_metadata = import_metadata(args.metadata)
    df_metadata = filter__metadata(df_metadata,df_counts)
    cond_dict = metadata__2__cond_dict(df_metadata,df_counts)

    print("Calculating gRNA fitness...")
    sample_set_list = fit_calc(df_counts, cond_dict, float(args.lower), float(args.upper), args.output)
    write_checkpoint(sample_set_list, args.output)


if __name__ == "__main__":
    main()