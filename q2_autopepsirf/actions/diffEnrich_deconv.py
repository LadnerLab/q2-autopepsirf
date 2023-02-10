from q2_pepsirf.format_types import (
    PepsirfLinkTSVFormat,
    PepsirfDMPFormat,
    PepsirfDeconvBatchDirFmt
)

import os

def diffEnrich_deconv(
        ctx,
        raw_data,
        bins,
        deconv_threshold,
        mapfile_suffix,
        outfile_suffix,
        linked,
        scoring_strategy="summation",
        infer_pairs_source=True,
        flexible_reps_source=False,
        s_enrich_source=False,
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
        score_filtering=False,
        score_tie_threshold=0.0,
        score_overlap_threshold=0.0,
        id_name_map=None,
        single_threaded=False,
        remove_file_types=False,
        pepsirf_binary="pepsirf"):

    diffEnrich = ctx.get_action("autopepsirf", "diffEnrich")
    deconv = ctx.get_action("pepsirf", "deconv_batch")

    (col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
     read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot, 
     zscore_scatter, colsum_scatter, zenrich_out
     ) = diffEnrich(
        raw_data=raw_data,
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

    (dir_out, score_per_round, map_dir) = deconv(
        enriched_dir=enrich_dir,
        threshold=deconv_threshold,
        mapfile_suffix=mapfile_suffix,
        outfile_suffix=outfile_suffix,
        linked=linked,
        scoring_strategy=scoring_strategy,
        score_filtering=score_filtering,
        score_tie_threshold=score_tie_threshold,
        score_overlap_threshold=score_overlap_threshold,
        id_name_map=id_name_map,
        single_threaded=single_threaded,
        remove_file_types=remove_file_types,
        outfile=os.path.join(pepsirf_tsv_dir, "deconv.out"),
        pepsirf_binary=pepsirf_binary
    )

    if pepsirf_tsv_dir and tsv_base_str:
        deconv_base = "%s_deconv_dir.tsv" % (tsv_base_str)
        deconv_tsv = dir_out.view(PepsirfDeconvBatchDirFmt)
        deconv_tsv.save(os.path.join(pepsirf_tsv_dir, deconv_base), ext=".tsv")    

    return (
        dir_out, score_per_round, map_dir, col_sum, diff, diff_ratio, 
        zscore_out, nan_out, sample_names, read_counts, rc_boxplot_out,
        enrich_dir, enrichedCountsBoxplot, zscore_scatter, colsum_scatter,
        zenrich_out
    )

