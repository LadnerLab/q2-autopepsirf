from q2_pepsirf.format_types import(
    PepsirfContingencyTSVFormat, PepsirfDemuxDiagnosticFormat
	)

import os

# Name: redoDemux
# Process: automatically runs through q2-ps-qc, q2-pepsirf, and q2_autopepsirf modules
# Method Input/Parameters:
# Method output/Returned:
# Dependencies:
# (autopepsirf: diffEnrich),
# (pepsirf: demux),
# (ps-qc: filter_counts_matrix, generate_corr_matrix)
def redoDemux(
	ctx,
	input_r1,
    index,
    samplelist,
    seq,
    bins,
    pairs_file=None,
    input_r2 = None,
    fif = None,
    library = None,
    read_per_loop = 100000,
    num_threads = 2,
    phred_base = 33,
    phred_min_score = 0,
    sindex = None,
    translate_aggregates = False,
    concatemer = False,
    sname = "SampleName",
    index1 = None,
    index2 = "0,0,0",
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
    correlation_threshold=0.8,
    pepsirf_binary="pepsirf",
    filtered_good_tsv_dir ="./filtered_good_tsv"
    ):

	if pepsirf_tsv_dir:
		if not os.path.isdir(pepsirf_tsv_dir):
			os.mkdir(pepsirf_tsv_dir)
	
	# collect actions
	demux = ctx.get_action("pepsirf", "demux")
	redoNoDemux = ctx.get_action("autopepsirf", "redoNoDemux")

	# run demux
	# TODO: make diagnostic data optional to speed up time
	(raw_counts, diagnostic_data) = demux(
							input_r1=input_r1,
				            index=index,
				            samplelist=samplelist,
				            seq=seq,
				            input_r2=input_r2,
				            fif=fif,
				            library=library,
				            read_per_loop=read_per_loop,
				            num_threads=num_threads,
				            phred_base=phred_base,
				            phred_min_score=phred_min_score,
				            sindex=sindex,
				            translate_aggregates=translate_aggregates,
				            concatemer=concatemer,
				            sname=sname,
				            index1=index1,
				            index2=index2,
				            outfile=os.path.join(pepsirf_tsv_dir, "demux.out"),
				            pepsirf_binary=pepsirf_binary
						)

	# run the rest with redoNoDemux
	(
        filtered_counts,
        bad_correlation_vis_zscores, good_correlation_vis_zscores,
		filtered_counts_good, col_sum_good, diff_good, diff_ratio_good, zscore_out_good,
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
    ) = redoNoDemux(
				    input_raw_counts=raw_counts,
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

	return (raw_counts, diagnostic_data,
        filtered_counts,
        bad_correlation_vis_zscores, good_correlation_vis_zscores,
		filtered_counts_good, col_sum_good, diff_good, diff_ratio_good, zscore_out_good,
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
    )