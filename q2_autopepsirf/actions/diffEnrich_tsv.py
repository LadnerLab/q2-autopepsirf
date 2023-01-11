import pandas as pd
import qiime2
from collections import defaultdict
import csv, os

from q2_pepsirf.format_types import (
    PepsirfInfoSumOfProbesFmt, PepsirfInfoSNPNFormat,
    PepsirfContingencyTSVFormat, ZscoreNanFormat,
    EnrichedPeptideDirFmt, PeptideBinFormat,
    EnrichThreshFileFormat
    )

# Name: diffenrich
# Process: automatically runs through q2-ps-plot modules and q2-pepsirf modules
# Method Input/Parameters: default ctx, raw_data_filepath, bins_filepath, negative_controls_filepath, negative_ids,
# negative_names, thresh_file_filepath, exact_z_thresh, exact_zenrich_thresh, step_z_thresh,
# upper_z_thresh, lower_z_thresh, raw_constraint, low_raw_reads, pepsirf_binary
# Method output/Returned: col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
# read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot, zscore_scatter, colsum_scatter
# Dependencies: (ps-plot: raedCountsBoxplot, enrichmentRCBoxplot, repScatters, zenrich), 
# (pepsirf: norm, zscore, infoSNPN, infoSumOfProbes, enrich)
def diffEnrich_tsv(
    ctx,
    raw_data_filepath,
    bins_filepath,
    infer_pairs_source=True,
    flexible_reps_source=False,
    s_enrich_source=False,
    user_defined_source = None,
    negative_control_filepath=None,
    negative_id=None,
    negative_names=None,
    thresh_file_filepath = None,
    exact_z_thresh = None,
    exact_cs_thresh = "20",
    exact_zenrich_thresh = None,
    pepsirf_tsv_dir = "./",
    tsv_base_str = None,
    step_z_thresh = 5,
    upper_z_thresh = 30,
    lower_z_thresh = 5,
    raw_constraint = 300000,
    low_raw_reads = False,
    hdi = 0.95,
    pepsirf_binary = "pepsirf"
    ):

    # collect diffEnrich action
    diffEnrich = ctx.get_action('autopepsirf', 'diffEnrich')

    # import raw data into an artifact
    raw_data = ctx.make_artifact(
        type='FeatureTable[RawCounts]',
        view=raw_data_filepath,
        view_type=PepsirfContingencyTSVFormat
    )

    # import bins into an artifact
    bins = ctx.make_artifact(
        type = 'PeptideBins',
        view = bins_filepath,
        view_type=PeptideBinFormat
    )

    # if negative_control provided import into artifact
    if negative_control_filepath:
        negative_control = ctx.make_artifact(
            type='FeatureTable[Normed]',
            view=negative_control_filepath,
            view_type=PepsirfContingencyTSVFormat
        )
    # otherwise set negative control to none
    else:
        negative_control = None
    
    #if thresh-file provided import into artifact
    if thresh_file_filepath:
        thresh_file = ctx.make_artifact(
            type='EnrichThresh',
            view=thresh_file_filepath,
            view_type=EnrichThreshFileFormat
        )
    #otherwise set thresh-file to none
    else:
        thresh_file = None

    # run the diffEnrich module with all the inputs/parameters given
    ( col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
    read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot, 
    zscore_scatter, colsum_scatter, zenrich_out, ) = diffEnrich(
        raw_data = raw_data,
        bins = bins,
        infer_pairs_source = infer_pairs_source,
        flexible_reps_source = flexible_reps_source,
        s_enrich_source = s_enrich_source,
        user_defined_source = user_defined_source,
        negative_control = negative_control,
        negative_id = negative_id,
        negative_names = negative_names,
        thresh_file = thresh_file,
        exact_z_thresh = exact_z_thresh,
        exact_cs_thresh = exact_cs_thresh,
        exact_zenrich_thresh = exact_zenrich_thresh,
        pepsirf_tsv_dir = pepsirf_tsv_dir,
        tsv_base_str = tsv_base_str,
        step_z_thresh = step_z_thresh,
        upper_z_thresh = upper_z_thresh,
        lower_z_thresh = lower_z_thresh,
        raw_constraint = raw_constraint,
        low_raw_reads = low_raw_reads,
        hdi = hdi,
        pepsirf_binary = pepsirf_binary 
    )

    # return all the qza and qzv files
    return (
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
        )
