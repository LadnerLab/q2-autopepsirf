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

def diffEnrich_tsv(
    ctx,
    raw_data_filepath,
    bins_filepath,
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
    hdi = 0.95,
    pepsirf_binary = "pepsirf"
    ):

    diffEnrich = ctx.get_action('autopepsirf', 'diffEnrich')

    raw_data = ctx.make_artifact(
        type='FeatureTable[RawCounts]',
        view=raw_data_filepath,
        view_type=PepsirfContingencyTSVFormat
    )

    bins = ctx.make_artifact(
        type = 'PeptideBins',
        view = bins_filepath,
        view_type=PeptideBinFormat
    )

    if negative_control_filepath:
        negative_control = ctx.make_artifact(
            type='FeatureTable[Normed]',
            view=negative_control_filepath,
            view_type=PepsirfContingencyTSVFormat
        )
    else:
        negative_control = None
    
    if thresh_file_filepath:
        thresh_file = ctx.make_artifact(
            type='EnrichThresh',
            view=thresh_file_filepath,
            view_type=EnrichThreshFileFormat
        )
    else:
        thresh_file = None

    ( col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
    read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot, 
    zscore_scatter, colsum_scatter, zenrich_out, ) = diffEnrich(
        raw_data = raw_data,
        bins = bins,
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
        hdi = hdi,
        pepsirf_binary = pepsirf_binary 
    )

    return (
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
        )