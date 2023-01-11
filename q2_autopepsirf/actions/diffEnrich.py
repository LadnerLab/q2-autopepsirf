from math import inf
import pandas as pd
import qiime2
from collections import defaultdict
import csv, os

from q2_pepsirf.format_types import PepsirfInfoSumOfProbesFmt, PepsirfInfoSNPNFormat, PepsirfContingencyTSVFormat, ZscoreNanFormat, EnrichedPeptideDirFmt

# Name: diffenrich
# Process: automatically runs through q2-ps-plot modules and q2-pepsirf modules
# Method Input/Parameters: default ctx, raw_data, bins, negative_controls, negative_ids,
# negative_names, thresh_file, exact_z_thresh, exact_zenrich_thresh, step_z_thresh,
# upper_z_thresh, lower_z_thresh, raw_constraint, low_raw_reads, pepsirf_binary
# Method output/Returned: col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
# read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot, zscore_scatter, colsum_scatter
# Dependencies: (ps-plot: raedCountsBoxplot, enrichmentRCBoxplot, repScatters, zenrich), 
# (pepsirf: norm, zscore, infoSNPN, infoSumOfProbes, enrich)
def diffEnrich(
    ctx,
    raw_data,
    bins,
    infer_pairs_source=True,
    flexible_reps_source=False,
    s_enrich_source=False,
    user_defined_source = None,
    negative_control=None,
    negative_id=None,
    negative_names=None,
    thresh_file = None,
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

    # if pepsirf_tsv_dir provided, make sure the provided dir is not a already created dir
    # otherwise, make it a dir
    if pepsirf_tsv_dir:
        if not os.path.isdir(pepsirf_tsv_dir):
            os.mkdir(pepsirf_tsv_dir)

        # if tsv_base_str not provided with pepsirf_tsv_dir set a default
        if not tsv_base_str:
            tsv_base_str = "aps-output"

    # collect the actions from ps-plot and q2-pepsirf to be executed
    norm = ctx.get_action('pepsirf', 'norm')
    zscore = ctx.get_action('pepsirf', 'zscore')
    infoSNPN = ctx.get_action('pepsirf', 'infoSNPN')
    enrich = ctx.get_action('pepsirf', 'enrich')
    infoSOP = ctx.get_action('pepsirf', 'infoSumOfProbes')
    RCBoxplot = ctx.get_action('ps-plot', 'readCountsBoxplot')
    enrichBoxplot = ctx.get_action('ps-plot', 'enrichmentRCBoxplot')
    repScatter = ctx.get_action('ps-plot', 'repScatters')
    zenrich = ctx.get_action('ps-plot', 'zenrich')

    # run norm module to recieved col-sum
    col_sum, = norm(peptide_scores = raw_data,
                    normalize_approach = "col_sum",
                    negative_control = None,
                    negative_id = None,
                    negative_names = None,
                    precision = 2,
                    outfile = os.path.join(pepsirf_tsv_dir, "norm.out"),
                    pepsirf_binary = pepsirf_binary)

    # convert the qza output into a tsv and save it
    if pepsirf_tsv_dir and tsv_base_str:
        cs_base = "%s_CS.tsv" % (tsv_base_str)
        cs_tsv = col_sum.view(PepsirfContingencyTSVFormat)
        cs_tsv.save(os.path.join(pepsirf_tsv_dir, cs_base), ext = ".tsv") #requires qiime2-2021.11

    # create list for collection of sample names
    if not negative_names and not negative_id:
        if not negative_control:
            negative_names = []
        else:
            negDF = negative_control.view(pd.DataFrame)
            negative_names = list(negDF.index)


    # run norm module to recieve diff
    diff, = norm(peptide_scores = col_sum,
                    normalize_approach = "diff",
                    negative_control = negative_control,
                    negative_id = negative_id,
                    negative_names = negative_names,
                    precision = 2,
                    outfile = os.path.join(pepsirf_tsv_dir, "norm.out"),
                    pepsirf_binary = pepsirf_binary)

    # convert the qza output into a tsv and save it
    if pepsirf_tsv_dir and tsv_base_str:
        diff_base = "%s_SBD.tsv" % (tsv_base_str)
        diff_tsv = diff.view(PepsirfContingencyTSVFormat)
        diff_tsv.save(os.path.join(pepsirf_tsv_dir, diff_base), ext = ".tsv")

    # run norm module to recieve diff-ratio
    diff_ratio, = norm(peptide_scores = col_sum,
                    normalize_approach = "diff_ratio",
                    negative_control = negative_control,
                    negative_id = negative_id,
                    negative_names = negative_names,
                    precision = 2,
                    outfile = os.path.join(pepsirf_tsv_dir, "norm.out"),
                    pepsirf_binary = pepsirf_binary)

    # convert the qza output into a tsv and save it
    if pepsirf_tsv_dir and tsv_base_str:
        diffR_base = "%s_SBDR.tsv" % (tsv_base_str)
        diffR_tsv = diff_ratio.view(PepsirfContingencyTSVFormat)
        diffR_tsv.save(os.path.join(pepsirf_tsv_dir, diffR_base), ext = ".tsv")

    # run zscore module to recieve zscore and nan files
    zscore_out, nan_out = zscore(
        scores = diff,
        bins = bins,
        hdi = hdi,
        outfile = os.path.join(pepsirf_tsv_dir, "zscore.out"),
        pepsirf_binary = pepsirf_binary
    )

    # convert the qza output into a tsv and save it
    if pepsirf_tsv_dir and tsv_base_str:
        zscore_base = "%s_Z-HDI%s.tsv" % (tsv_base_str, str(int(hdi*100)))
        zscore_tsv = zscore_out.view(PepsirfContingencyTSVFormat)
        zscore_tsv.save(os.path.join(pepsirf_tsv_dir, zscore_base), ext = ".tsv")

        nan_base = "%s_Z-HDI%s.nan" % (tsv_base_str, str(int(hdi*100)))
        nan_tsv = nan_out.view(ZscoreNanFormat)
        nan_tsv.save(os.path.join(pepsirf_tsv_dir, nan_base), ext = ".nan")


    # run info module to collect sample names
    sample_names, = infoSNPN(
        input = raw_data,
        get = "samples",
        outfile = os.path.join(pepsirf_tsv_dir, "info.out"),
        pepsirf_binary = pepsirf_binary
    )

    # convert the qza output into a tsv and save it
    if pepsirf_tsv_dir and tsv_base_str:
        sn_base = "%s_SN.tsv" % (tsv_base_str)
        sn_tsv = sample_names.view(PepsirfInfoSNPNFormat)
        sn_tsv.save(os.path.join(pepsirf_tsv_dir, sn_base), ext = ".tsv")

    # run info to collect read counts
    read_counts, = infoSOP(
        input = raw_data,
        outfile = os.path.join(pepsirf_tsv_dir, "info.out"),
        pepsirf_binary = pepsirf_binary
    )

    # convert the qza output into a tsv and save it
    if pepsirf_tsv_dir and tsv_base_str:
        rc_base = "%s_RC.tsv" % (tsv_base_str)
        rc_tsv = read_counts.view(PepsirfInfoSumOfProbesFmt)
        rc_tsv.save(os.path.join(pepsirf_tsv_dir, rc_base), ext = ".tsv")

    # run readCounts boxplot module to recieve visualization
    rc_boxplot_out, = RCBoxplot(
        read_counts = read_counts,
        png_out_dir = pepsirf_tsv_dir
    )

    # create variables for source file creation
    if infer_pairs_source or flexible_reps_source or s_enrich_source:
        sourceDic = defaultdict(list)
        sampleNM = sample_names.view(PepsirfInfoSNPNFormat)
        source = os.path.join(pepsirf_tsv_dir, "samples_source.tsv")


        # open samples file and collect samples into a dictionary
        with open( str(sampleNM) ) as SN:
            for line in SN:
                sample = line.strip()
                sourceLS = sample.rsplit('_', 1)
                sourced = sourceLS[0]
                sourceDic[sourced].append(sample)
                if not negative_names and not negative_id and not negative_control:
                    negative_names.append(sample)

        # create a source file written with column 1 as the sample names
        # and the column 2 as the source column
        # the source file will be put in the tsv directory
        with open( source , "w" ) as tsvWriter:
            writer = csv.writer(tsvWriter, delimiter='\t')
            writer.writerow(['sampleID', 'source'])
            for srce, samples in sourceDic.items():
                if flexible_reps_source:
                    for name in samples:
                            writer.writerow([name, srce])
                elif s_enrich_source:
                    for name in samples:
                            writer.writerow([name, name])
                elif infer_pairs_source:
                    if len(samples) > 1:
                        for name in samples:
                            writer.writerow([name, srce])

        # convert source file to metadata column to be used within the modules
        source_col = qiime2.Metadata.load(source).get_column("source")
    
    elif user_defined_source:
        source_col = user_defined_source

    # run enrich module
    enrich_dir, = enrich(
        source = source_col,
        thresh_file = thresh_file,
        zscores = zscore_out,
        col_sum = col_sum,
        exact_z_thresh = exact_z_thresh,
        exact_cs_thresh = exact_cs_thresh,
        raw_scores = raw_data,
        raw_constraint = raw_constraint,
        enrichment_failure = True,
        low_raw_reads = low_raw_reads,
        outfile = os.path.join(pepsirf_tsv_dir, "enrich.out"),
        pepsirf_binary = pepsirf_binary
    )

    # convert the qza output into a tsv and save it
    if pepsirf_tsv_dir and tsv_base_str:
        if exact_z_thresh:
            enrich_zt = exact_z_thresh.split(",")
            enrich_cst = exact_cs_thresh.split(",")
            if len(enrich_zt) > 1:
                enrich_base = "%s-%sZ-HDI%s_" % (
                    enrich_zt[0], enrich_zt[1], str(int(hdi*100)))
            else:
                enrich_base = "%sZ-HDI%s_" % (
                    enrich_zt[0], str(int(hdi*100)))
            if len(enrich_cst) > 1:
                enrich_base += "%s-%sCS_%sraw" % (
                    enrich_cst[0], enrich_cst[1], str(raw_constraint)
                )
            else:
                enrich_base += "%sCS_%sraw" % (
                    enrich_cst[0], str(raw_constraint)
                )
        else:
            enrich_base = "enriched"
        enrich_tsv = enrich_dir.view(EnrichedPeptideDirFmt)
        enrich_tsv.save(os.path.join(pepsirf_tsv_dir, enrich_base))

    # run enrichment boxplot module to recieve visualization
    enrichedCountsBoxplot, = enrichBoxplot(
        enriched_dir = enrich_dir,
        png_out_dir = pepsirf_tsv_dir
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

    # run the zenrich module to collect visualization
    zenrich_out, = zenrich(
        data = col_sum,
        zscores = zscore_out,
        negative_controls = negative_names,
        negative_id = negative_id,
        source = source_col,
        negative_data = negative_control,
        step_z_thresh = step_z_thresh,
        upper_z_thresh = upper_z_thresh,
        lower_z_thresh = lower_z_thresh,
        exact_z_thresh = exact_zenrich_thresh,
        exact_cs_thresh = exact_cs_thresh,
        pepsirf_binary = pepsirf_binary
    )

    # return all files created
    return (
        col_sum, diff, diff_ratio, zscore_out, nan_out, sample_names,
        read_counts, rc_boxplot_out, enrich_dir, enrichedCountsBoxplot,
        zscore_scatter, colsum_scatter, zenrich_out
        )
