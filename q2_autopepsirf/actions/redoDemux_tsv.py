from q2_pepsirf.format_types import(
    PepsirfDemuxFastqFmt, PepsirfDemuxIndexFmt,
    PepsirfDemuxSampleListFmt, PepsirfDemuxFifFmt,
    PepsirfDemuxLibraryFmt,
    PepsirfInfoSumOfProbesFmt, PepsirfInfoSNPNFormat,
    PepsirfContingencyTSVFormat, ZscoreNanFormat, EnrichedPeptideDirFmt,
    PeptideBinFormat, EnrichThreshFileFormat
	)

def redoDemux_tsv(
    ctx,
    input_r1_filepath,
    index_filepath,
    samplelist_filepath,
    seq,
    bins_filepath,
    pairs_file=None,
    input_r2_filepath = None,
    fif_filepath = None,
    library_filepath= None,
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
    infer_pairs_source=True,
    flexible_reps_source=False,
    s_enrich_source=False,
    user_defined_source=None,
    negative_control_filepath=None,
    negative_id=None,
    negative_names=None,
    thresh_file_filepath=None,
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
    
    # collect redoDemux action
    redoDemux = ctx.get_action("autopepsirf", "redoDemux")

    # import required files into artifacts
    input_r1 = ctx.make_artifact(
            type="DemuxFastq",
            view=input_r1_filepath,
            view_type=PepsirfDemuxFastqFmt
        )

    index = ctx.make_artifact(
            type="DemuxIndex",
            view=index_filepath,
            view_type=PepsirfDemuxIndexFmt
        )

    samplelist = ctx.make_artifact(
            type="DemuxSampleList",
            view=samplelist_filepath,
            view_type=PepsirfDemuxSampleListFmt
        )
    bins = ctx.make_artifact(
        type="PeptideBins",
        view=bins_filepath,
        view_type=PeptideBinFormat
    )


    # import optional files into artifacts
    if input_r2_filepath:
        input_r2 = ctx.make_artifact(
            type="DemuxFastq",
            view=input_r2_filepath,
            view_type=PepsirfDemuxFastqFmt
        )
    else:
        input_r2 = None

    if fif_filepath:
        fif = ctx.make_artifact(
            type="DemuxFif",
            view=fif_filepath,
            view_type=PepsirfDemuxFifFmt
        )
    else:
        fif = None

    if library_filepath:
        library = ctx.make_artifact(
            type="DemuxLibrary",
            view=library_filepath,
            view_type=PepsirfDemuxLibraryFmt
        )
    else:
        library = None

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

    # run redoDemux
    (raw_counts, diagnostic_data,
        filtered_counts, 
        bad_correlation_vis_col_sum, good_correlation_vis_col_sum,
        bad_correlation_vis_zscores, good_correlation_vis_zscores,
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
    ) = redoDemux(
            input_r1=input_r1,
            index=index,
            samplelist=samplelist,
            seq=seq,
            bins=bins,
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
            outfile=outfile,
            count_thresh=count_thresh,
            max_zeros=max_zeros,
            drop_samp_out=drop_samp_out,
            pairs_file=pairs_file,
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
            log_normalization=log_normalization,
            correlation_threshold=correlation_threshold,
            pepsirf_binary="pepsirf"
        )

    return (raw_counts, diagnostic_data,
        filtered_counts, 
        bad_correlation_vis_col_sum, good_correlation_vis_col_sum,
        bad_correlation_vis_zscores, good_correlation_vis_zscores,
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
    )






