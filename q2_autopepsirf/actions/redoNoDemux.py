from q2_pepsirf.format_types import PepsirfContingencyTSVFormat
import tempfile
import os

def redoNoDemux(
    ctx,
    input_raw_counts,
    bins,
    pairs_file=None,
    count_thresh=None,
    max_zeros=None,
    drop_samp_out="filtered_out.tsv",
    infer_pairs_source=False,
    flexible_reps_source=False,
    s_enrich_source=True,
    user_defined_source = None,
    negative_control=None,
    negative_id=None,
    negative_names=None,
    thresh_file=None,
    exact_z_thresh=None,
    exact_cs_thresh="20",
    exact_zenrich_thresh=None,
    pepsirf_tsv_dir="./",
    tsv_base_str=None,
    step_z_thresh=5,
    upper_z_thresh=30,
    lower_z_thresh=5,
    raw_constraint=300000,
    hdi=0.95,
    log_normalization=False,
    correlation_threshold=0.8,
    pepsirf_binary="pepsirf"
    ):
    
    filter_counts_matrix = ctx.get_action("q2-ps-qc", "filter_counts_matrix")
    diffEnrich = ctx.get_action("autopepsirf", "diffEnrich")
    generate_corr_matrix = ctx.get_action("q2-ps-qc", "generate_corr_matrix")

    # filter raw counts
    filtered_counts, = filter_counts_matrix(
                    input_matrix=input_raw_counts,
                    count_thresh=count_thresh,
                    max_zeros=max_zeros,
                    drop_samp_out=drop_samp_out
                    )

    with tempfile.TemporaryDirectory() as temp_dir:
        filtered_counts_filepath = os.path.join(temp_dir, "filtered_counts_matrix.tsv")
        filtered_counts.view(PepsirfContingencyTSVFormat).save(filtered_counts_filepath, ext=".tsv")

        # generated correlated matrix for col sum normalized matrix
        (bad_correlation_vis_filtered_counts, good_correlation_vis_filtered_counts
            ) = generate_corr_matrix(
                data=filtered_counts_filepath,
                samples=pairs_file,
                log_normalization=log_normalization,
                correlation_threshold=correlation_threshold,
                bad_corr_out="bad_corr_filtered_count.tsv",
                good_corr_out="good_corr_filtered_counts.tsv"
        )

        good_corr_filtered_counts = ctx.make_artifact(
            type="FeatureTable[RawCounts]",
            view="good_corr_filtered_counts.tsv",
            view_type=PepsirfContingencyTSVFormat
        )

    # run diffEnrich for zscores and normalized counts
    (col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
     read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot, 
     zscore_scatter, colsum_scatter, zenrich_out
     ) = diffEnrich(
        raw_data=good_corr_filtered_counts,
        bins=bins,
        infer_pairs_source=infer_pairs_source,
        flexible_reps_source=flexible_reps_source,
        s_enrich_source=s_enrich_source,
        user_defined_source=user_defined_source,
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
        pepsirf_binary=pepsirf_binary 
    )

    return (
        filtered_counts, 
        bad_correlation_vis_filtered_counts, good_correlation_vis_filtered_counts,
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
    )