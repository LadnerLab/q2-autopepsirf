from q2_pepsirf.format_types import (
    PepsirfContingencyTSVFormat, PepsirfInfoSumOfProbesFmt, PepsirfInfoSNPNFormat,
    PepsirfContingencyTSVFormat, ZscoreNanFormat, EnrichedPeptideDirFmt,
    PeptideBinFormat, EnrichThreshFileFormat
)

def redoNoDemux_tsv(
    ctx,
    input_raw_counts_filepath,
    bins_filepath,
    pairs_file=None,
    count_thresh=None,
    max_zeros=None,
    drop_samp_out="filtered_out.tsv",
    infer_pairs_source=False,
    flexible_reps_source=False,
    s_enrich_source=True,
    user_defined_source=None,
    negative_control_filepath=None,
    negative_id=None,
    negative_names=None,
    thresh_file_filepath=None,
    exact_z_thresh=None,
    exact_cs_thresh="20",
    exact_zenrich_thresh=None,
    pepsirf_tsv_dir="./pepsirf_tsv",
    tsv_base_str=None,
    step_z_thresh=5,
    upper_z_thresh=30,
    lower_z_thresh=5,
    raw_constraint=300000,
    hdi=0.95,
    log_normalization=False,
    correlation_threshold=0.8,
    pepsirf_binary="pepsirf",
    filtered_good_tsv_dir ="./filtered_good_tsv"
    ):
    
    redoNoDemux = ctx.get_action("autopepsirf", "redoNoDemux")

    # create artifacts
    input_raw_counts = ctx.make_artifact(
        type="FeatureTable[RawCounts]",
        view=input_raw_counts_filepath,
        view_type=PepsirfContingencyTSVFormat
    )

    bins = ctx.make_artifact(
        type="PeptideBins",
        view=bins_filepath,
        view_type=PeptideBinFormat
    )

    # create optional artifacts
    # if negative_control provided import into artifact
    if negative_control_filepath:
        negative_control = ctx.make_artifact(
            type="FeatureTable[Normed]",
            view=negative_control_filepath,
            view_type=PepsirfContingencyTSVFormat
        )
    # otherwise set negative control to none
    else:
        negative_control = None
    
    #if thresh-file provided import into artifact
    if thresh_file_filepath:
        thresh_file = ctx.make_artifact(
            type="EnrichThresh",
            view=thresh_file_filepath,
            view_type=EnrichThreshFileFormat
        )
    #otherwise set thresh-file to none
    else:
        thresh_file = None

    (
        filtered_counts,
        bad_correlation_vis_zscores, good_correlation_vis_zscores,
        filtered_counts_good, col_sum_good, diff_good, diff_ratio_good, zscore_out_good,
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
    ) = redoNoDemux(
            input_raw_counts=input_raw_counts,
            bins=bins,
            count_thresh=count_thresh,
            max_zeros=max_zeros,
            drop_samp_out=drop_samp_out,
            pairs_file=pairs_file,
            infer_pairs_source=infer_pairs_source,
            flexible_reps_source=flexible_reps_source,
            s_enrich_source=s_enrich_source,
            user_defined_source = user_defined_source,
            negative_control=negative_control,
            negative_id=negative_id,
            negative_names=negative_names,
            thresh_file=thresh_file,
            exact_z_thresh=exact_z_thresh,
            exact_cs_thresh=exact_cs_thresh,
            exact_zenrich_thresh=exact_zenrich_thresh,
            pepsirf_tsv_dir=pepsirf_tsv_dir,
            tsv_base_str=tsv_base_str,
            step_z_thresh=step_z_thresh,
            upper_z_thresh=upper_z_thresh,
            lower_z_thresh=lower_z_thresh,
            raw_constraint=raw_constraint,
            hdi=hdi,
            log_normalization=log_normalization,
            correlation_threshold=correlation_threshold,
            pepsirf_binary=pepsirf_binary,
            filtered_good_tsv_dir=filtered_good_tsv_dir
        )

    return (
        filtered_counts,
        bad_correlation_vis_zscores, good_correlation_vis_zscores,
        filtered_counts_good, col_sum_good, diff_good, diff_ratio_good, zscore_out_good,
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
    )