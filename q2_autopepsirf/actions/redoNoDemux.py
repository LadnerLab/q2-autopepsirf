from q2_pepsirf.format_types import PepsirfContingencyTSVFormat
import tempfile
import os
import pandas as pd

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
    pepsirf_tsv_dir="./pepsirf_tsv",
    tsv_base_str=None,
    step_z_thresh=5,
    upper_z_thresh=30,
    lower_z_thresh=5,
    raw_constraint=300000,
    hdi=0.95,
    log_normalization=False,
    reused_samples_in_pairs=True,
    correlation_threshold=0.8,
    pepsirf_binary="pepsirf",
    filtered_good_tsv_dir ="./filtered_good_tsv"
    ):
    
    if filtered_good_tsv_dir:
        if not os.path.isdir(filtered_good_tsv_dir):
            os.mkdir(filtered_good_tsv_dir)

    filter_counts_matrix = ctx.get_action("q2-ps-qc", "filter_counts_matrix")
    diffEnrich = ctx.get_action("autopepsirf", "diffEnrich")
    generate_corr_matrix = ctx.get_action("q2-ps-qc", "generate_corr_matrix")

    # filter raw counts based on counts and max zeros
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

    # get good correlated zscores
    with tempfile.TemporaryDirectory() as temp_dir:
        zscore_matrix_filepath = os.path.join(temp_dir, "zscore_matrix.tsv")
        zscore_out.view(PepsirfContingencyTSVFormat).save(zscore_matrix_filepath, ext=".tsv")

        # generated correlated matrix for zscores
        (bad_correlation_vis_zscores, good_correlation_vis_zscores
            ) = generate_corr_matrix(
                data=zscore_matrix_filepath,
                samples=pairs_file,
                log_normalization=log_normalization,
                reused_samples_in_pairs=reused_samples_in_pairs,
                correlation_threshold=correlation_threshold,
                bad_corr_out="bad_corr_zscores.tsv",
                good_corr_out="good_corr_zscores.tsv",
                bad_pairs_out="bad_pairs.tsv",
                good_pairs_out="good_pairs.tsv"
        )
    
        # filter raw counts based on good zscore correlatons
        good_samples = pd.read_csv("good_corr_zscores.tsv", sep="\t", index_col=0).columns.to_list()

    filtered_counts_good = filter_matrix(ctx, filtered_counts, good_samples, "FeatureTable[RawCounts]", filtered_good_tsv_dir, "filtered_counts_good.tsv")
    col_sum_good = filter_matrix(ctx, col_sum, good_samples, "FeatureTable[Normed]", filtered_good_tsv_dir, "col_sum_good.tsv")
    diff_good = filter_matrix(ctx, diff, good_samples, "FeatureTable[NormedDifference]", filtered_good_tsv_dir, "diff_good.tsv")
    diff_ratio_good = filter_matrix(ctx, diff_ratio, good_samples, "FeatureTable[NormedDiffRatio]", filtered_good_tsv_dir, "diff_ratio_good.tsv")
    zscore_out_good = filter_matrix(ctx, zscore_out, good_samples, "FeatureTable[Zscore]", filtered_good_tsv_dir, "zscore_out_good.tsv")
        


    return (
        filtered_counts,
        bad_correlation_vis_zscores, good_correlation_vis_zscores,
        filtered_counts_good, col_sum_good, diff_good, diff_ratio_good, zscore_out_good,
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
    )


def filter_matrix(ctx, matrix_artifact, good_samples, type, out_dir, filename):
    matrix_filepath = os.path.join(out_dir, filename)
    matrix_artifact.view(PepsirfContingencyTSVFormat).save(matrix_filepath, ext=".tsv")
    matrix_df = pd.read_csv(matrix_filepath, sep="\t", index_col=0)
    filtered_matrix_df = matrix_df[good_samples]
    filtered_matrix_df.to_csv(matrix_filepath, sep="\t")
    filtered_matrix_artifact = ctx.make_artifact(
            type=type,
            view=matrix_filepath,
            view_type=PepsirfContingencyTSVFormat
        )

    return filtered_matrix_artifact
