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

    # run diffEnrich for zscores and normalized counts
    (col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
     read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot, 
     zscore_scatter, colsum_scatter, zenrich_out
     ) = diffEnrich(
        raw_data=filtered_counts,
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

    with tempfile.TemporaryDirectory() as temp_dir:
        zscore_matrix_filepath = os.path.join(temp_dir, "zscore_matrix.tsv")
        zscore_out.view(PepsirfContingencyTSVFormat).save(zscore_matrix_filepath, ext=".tsv")

        col_sum_matrix_filepath = os.path.join(temp_dir, "col_sum_matrix.tsv")
        col_sum.view(PepsirfContingencyTSVFormat).save(col_sum_matrix_filepath, ext=".tsv")

        
        # generated correlated matrix for col sum normalized matrix
        (bad_correlation_vis_col_sum, good_correlation_vis_col_sum
            ) = generate_corr_matrix(
                data=col_sum_matrix_filepath,
                samples=pairs_file,
                log_normalization=log_normalization,
                correlation_threshold=correlation_threshold,
                bad_corr_out="bad_corr_col_sum.tsv",
                good_corr_out="good_corr_col_sum.tsv"
        )

        # generated correlated matrix for zscores
        (bad_correlation_vis_zscores, good_correlation_vis_zscores
            ) = generate_corr_matrix(
                data=zscore_matrix_filepath,
                samples=pairs_file,
                log_normalization=log_normalization,
                correlation_threshold=correlation_threshold,
                bad_corr_out="bad_corr_zscores.tsv",
                good_corr_out="good_corr_zscores.tsv"
        )

    return (
        filtered_counts, 
        bad_correlation_vis_col_sum, good_correlation_vis_col_sum,
        bad_correlation_vis_zscores, good_correlation_vis_zscores,
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
    )