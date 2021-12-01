import pandas as pd
import qiime2
from collections import defaultdict
import csv, os

from q2_pepsirf.format_types import PepsirfInfoSNPNFormat

# Name: diffenrich
# Process: automatically runs through q2-ps-plot modules and q2-pepsirf modules
# Method Input/Parameters: default ctx, raw_data, bins, negative_controls, negative_ids,
# negative_names, thresh_file, exact_z_thresh, raw_constraint, pepsirf_binary
# Method output/Returned: col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
# read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot, zscore_scatter, colsum_scatter
# Dependencies: (ps-plot: raedCountsBoxplot, enrichmentRCBoxplot, repScatters), 
# (pepsirf: norm, zscore, infoSNPN, infoSumOfProbes, enrich)
def diffEnrich(
    ctx,
    raw_data,
    bins,
    negative_control=None,
    negative_id=None,
    negative_names=None,
    thresh_file = None,
    exact_z_thresh = None,
    raw_constraint = 300000,
    pepsirf_binary = "pepsirf"
):

    # collect the actions from ps-plot and q2-pepsirf to be executed
    norm = ctx.get_action('pepsirf', 'norm')
    zscore = ctx.get_action('pepsirf', 'zscore')
    infoSNPN = ctx.get_action('pepsirf', 'infoSNPN')
    enrich = ctx.get_action('pepsirf', 'enrich')
    infoSOP = ctx.get_action('pepsirf', 'infoSumOfProbes')
    RCBoxplot = ctx.get_action('ps-plot', 'readCountsBoxplot')
    enrichBoxplot = ctx.get_action('ps-plot', 'enrichmentRCBoxplot')
    repScatter = ctx.get_action('ps-plot', 'repScatters')

    # run norm module to recieved col-sum
    col_sum, = norm(peptide_scores = raw_data,
                    normalize_approach = "col_sum",
                    negative_control = None,
                    negative_id = None,
                    negative_names = None,
                    precision = 2,
                    pepsirf_binary = pepsirf_binary)

    # run norm module to recieve diff
    diff, = norm(peptide_scores = col_sum,
                    normalize_approach = "diff",
                    negative_control = negative_control,
                    negative_id = negative_id,
                    negative_names = negative_names,
                    precision = 2,
                    pepsirf_binary = pepsirf_binary)

    # run norm module to recieve diff-ratio
    diff_ratio, = norm(peptide_scores = col_sum,
                    normalize_approach = "diff_ratio",
                    negative_control = negative_control,
                    negative_id = negative_id,
                    negative_names = negative_names,
                    precision = 2,
                    pepsirf_binary = pepsirf_binary)

    # run zscore module to recieve zscore and nan files
    zscore_out, nan_out = zscore(
        scores = diff,
        bins = bins,
        hdi = 0.95,
        pepsirf_binary = pepsirf_binary
    )

    # run info module to collect sample names
    sample_names, = infoSNPN(
        input = raw_data,
        get = "samples",
        pepsirf_binary = pepsirf_binary
    )

    # run info to collect read counts
    read_counts, = infoSOP(
        input = raw_data,
        pepsirf_binary = pepsirf_binary
    )

    # run readCounts boxplot module to recieve visualization
    rc_boxplot_out, = RCBoxplot(
        read_counts = read_counts
    )

    # create variables for source file creation
    sourceDic = defaultdict(list)
    sampleNM = sample_names.view(PepsirfInfoSNPNFormat)
    source = "pairs_source.tsv"

    # open samples file and collect samples into a dictionary
    with open( str(sampleNM) ) as SN:
        for line in SN:
            sample = line.strip()
            sourceLS = sample.rsplit('_', 1)
            sourced = sourceLS[0]
            sourceDic[sourced].append(sample)

    # create a source file written with column 1 as the sample names
    # and the column 2 as the source column
    with open( source , "w" ) as tsvWriter:
        writer = csv.writer(tsvWriter, delimiter='\t')
        writer.writerow(['sampleID', 'source'])
        for srce, samples in sourceDic.items():
            for name in samples:
                writer.writerow([name, srce])

    # convert source file to metadata column to be used within the modules
    source_col = qiime2.Metadata.load(source).get_column("source")

    # run enrich module
    enrich_dir, = enrich(
        source = source_col,
        thresh_file = thresh_file,
        zscores = zscore_out,
        exact_z_thresh = exact_z_thresh,
        raw_scores = raw_data,
        raw_constraint = raw_constraint,
        enrichment_failure = True,
        pepsirf_binary = pepsirf_binary
    )

    # run enrichment boxplot module to recieve visualization
    enrichedCountsBoxplot, = enrichBoxplot(
        enriched_dir = enrich_dir
    )

    # run repScatter module to collect visualization
    zscore_scatter, = repScatter(
        source = source_col,
        plot_log = False,
        zscore = zscore_out
    )

    # run repScatter module to collect visualization
    colsum_scatter, = repScatter(
        source = source_col,
        plot_log = True,
        col_sum = col_sum
    )

    # return all files created
    return (
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter
        )