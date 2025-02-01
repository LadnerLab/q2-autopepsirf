from q2_autopepsirf.actions.diffEnrich import diffEnrich
from q2_autopepsirf.actions.diffEnrich_tsv import diffEnrich_tsv
from q2_autopepsirf.actions.diffEnrich_deconv import diffEnrich_deconv
from q2_autopepsirf.actions.diffEnrich_deconv_tsv import diffEnrich_deconv_tsv
from q2_autopepsirf.actions.redoDemux import redoDemux
from q2_autopepsirf.actions.redoDemux_tsv import redoDemux_tsv
from q2_autopepsirf.actions.redoNoDemux import redoNoDemux
from q2_autopepsirf.actions.redoNoDemux_tsv import redoNoDemux_tsv

from q2_types.feature_table import FeatureTable
from qiime2.plugin import (
    Plugin, TypeMap, Str, List, MetadataColumn,
    Categorical, Int, Range, Visualization, Float,
    Bool
)
from q2_pepsirf.format_types import (
    RawCounts, Normed, NormedDifference,
    NormedDiffRatio, PeptideBins, Zscore,
    ZscoreNan, InfoSNPN, EnrichThresh,
    PairwiseEnrichment, InfoSumOfProbes,
    DeconvBatch, PeptideAssignmentMap,
    ScorePerRound, Link, PepsirfDMP,
    DemuxFastq, DemuxIndex, DemuxSampleList,
    DemuxFif, DemuxLibrary , DemuxDiagnostic
)

import importlib
import q2_autopepsirf

# This is the plugin object. It is what the framework will load and what an
# interface will interact with. Basically every registration we perform will
# involve this object in some way.
plugin = Plugin(
    "autopepsirf",
    version=q2_autopepsirf.__version__,
    website="https://github.com/LadnerLab/q2-autopepsirf",
    description="Qiime2 plugin used for the automation of q2-pepsirf and q2-ps-plot."
)

# shared outputs for diffEnrich and diffEnrich tsv pipeline
diffEnrich_shared_outputs = [
    ("col_sum", FeatureTable[Normed]),
    ("diff", FeatureTable[NormedDifference]),
    ("diff_ratio", FeatureTable[NormedDiffRatio]),
    ("zscore", FeatureTable[Zscore]),
    ("zscore_nan", ZscoreNan),
    ("sample_names", InfoSNPN),
    ("read_counts", InfoSumOfProbes),
    ("rc_boxplot", Visualization),
    ("enrich", PairwiseEnrichment),
    ("enrich_count_boxplot", Visualization),
    ("zscore_scatter", Visualization),
    ("colsum_scatter", Visualization),
    ("zenrich_scatter", Visualization)
]

diffEnrich_good_outputs = [
    ("filtered_counts_good", FeatureTable[RawCounts]),
    ("col_sum_good", FeatureTable[Normed]),
    ("diff_good", FeatureTable[NormedDifference]),
    ("diff_ratio_good", FeatureTable[NormedDiffRatio]),
    ("zscore_good", FeatureTable[Zscore])
]

# shared paremters for diffEnrich and diffEnrich tsv pipeline
diffEnrich_shared_parameters = {
    "negative_id": Str,
    "negative_names": List[Str],
    "pepsirf_binary": Str,
    "exact_z_thresh": Str,
    "exact_cs_thresh": Str,
    "raw_constraint": Int % Range(0, None),
    "exact_zenrich_thresh": List[Str],
    "step_z_thresh": Int % Range(1, None),
    "upper_z_thresh": Int % Range(2, None),
    "lower_z_thresh": Int % Range(1, None),
    "pepsirf_tsv_dir": Str,
    "tsv_base_str": Str,
    "hdi": Float % Range(0.0, 1.0),
    "infer_pairs_source": Bool,
    "flexible_reps_source": Bool,
    "s_enrich_source": Bool,
    "user_defined_source": MetadataColumn[Categorical]
}

# shared parameter descriptions for diffEnrich and diffEnrich tsv pipeline
diffEnrich_shared_parameter_description = {
    "negative_id": "Optional approach for identifying negative controls."
        " Provide a unique string at the start of all negative control"
        " samples.",
    "negative_names": "Optional approach for identifying negative controls."
        " Space-separated list of negative control sample names.",
    "pepsirf_binary": "The binary to call pepsirf on your system.",
    "exact_z_thresh": "Individual Exact z score threshold separated by a comma"
        " for creation of threshold file to run pepsirf's enrich module"
        " (Ex: 6,10 or 30)",
    "exact_cs_thresh": "Individual Exact col-sum threshold separated by a"
        " comma for creation of threshold file to run pepsirf's enrich module"
        " (Ex: 6,10 or 30)",
    "raw_constraint": "The minimum total raw count across all peptides for a"
        " sample to be included in the analysis. This provides a way to impose"
        " a minimum read count for a sample to be evaluated.",
    "exact_zenrich_thresh": "List of exact z score thresholds either"
        " individual or combined. List MUST BE in descending order. (Example"
        " argument: '--p-exact-zenrich-thresh 25 10 3' or"
        " '--p-exact-zenrich-thresh 6,25 4,10 1,3')",
    "step_z_thresh": "Integar to increment z-score thresholds.",
    "upper_z_thresh": "Upper limit of z-score thresholds (non-inclusive).",
    "lower_z_thresh": "Lower limit of z-score thresholds (inclusive).",
    "pepsirf_tsv_dir": "Provide a directory path. Must also provide"
        " tsv-base-str for output of tsv verison of qza files. The"
        " source_samples file and png boxplot outputs will always be put"
        " within this directory.",
    "tsv_base_str": "The base name for the output tsv files excluding ay"
        " extensions, typcally the raw data filename (EX: --p-tsv-base-str"
        " raw_data). Must also provide pepsirf-tsv-dir, if pepsirf-tsv-dir"
        " provided without tsv-base-str, the default will be 'aps-output'.",
    "hdi": "Alternative approach for discarding outliers prior to calculating"
        " mean and stdev. If provided, this argument will override --trim,"
        " which trims evenly from both sides of the distribution. For --hdi,"
        " the user should provide the high density interval to be used for"
        " calculation of mean and stdev. For example, '--hdi 0.95' would"
        " instruct the program to utilize the 95% highest density interval"
        " (from each bin) for these calculations.",
    "infer_pairs_source": "Infer sample pairs from names. This option assumes"
        " names of replicates will be identical with the exception of a final"
        " string denoted with a '_'. For example, these names would be"
        " considered two replicates of the same sample: VW_100_1X_A and"
        " VW_100_1X_B",
    "flexible_reps_source": "Will infer the number of replicates for each"
        " sample based on sample names, and will not require any specific"
        " number of replicates for inclusion. Therefore, some samples may have"
        " a single replicate, some may have 2, 3, 4 etc. And all replicates of"
        " a given sample will be considered for determining enriched"
        " peptides.",
    "s_enrich_source": "All samples will be processed individually as samples"
        " with only one replicate",
    "user_defined_source": "Metadata file containing all sample names and"
        " their source groups. Used to create pairs tsv to run pepsirf enrich"
        " module."
}

diffEnrich_inputs = {
    "negative_control": FeatureTable[Normed],
    "bins": PeptideBins,
    "thresh_file": EnrichThresh
    }

diffEnrich_input_descriptions = {
    "negative_control": "Name of FeatureTable matrix file containing data"
        " for sb samples.",
    "bins": "Name of the file containing bins, one bin per line, as output"
        " by the bin module. Each bin contains a tab-delimited list of"
        " peptide names.",
    "thresh_file": "The name of a tab-delimited file containing one"
        " tab-delimited matrix filename and threshold(s), one per line. If"
        " providing more than z score matrix."
}

diffEnrich_tsv_params = {
    "negative_control_filepath": Str,
    "bins_filepath": Str,
    "thresh_file_filepath": Str,
    }

diffEnrich_tsv_param_descriptions = {
    "negative_control_filepath": "Name of .tsv matrix file containing data"
        " for sb samples.",
    "bins_filepath": "Name of the file containing bins, one bin per line,"
        " as output by the bin module. Each bin contains a tab-delimited"
        " list of peptide names.",
    "thresh_file_filepath": "The name of a tab-delimited file containing"
        " one tab-delimited matrix filename and threshold(s), one per"
        " line. If providing more than z score matrix."
}

# action set up for diffEnrich module
plugin.pipelines.register_function(
    function=diffEnrich,
    inputs={
        "raw_data": FeatureTable[RawCounts],
        **diffEnrich_inputs
    },
    outputs=diffEnrich_shared_outputs,
    parameters=diffEnrich_shared_parameters,
    input_descriptions={
        "raw_data": "Raw data matrix.",
        **diffEnrich_input_descriptions
    },
    output_descriptions=None,
    parameter_descriptions=diffEnrich_shared_parameter_description,
    name="diffEnrich Pepsirf Pipeline",
    description="Uses the diff normaization from pepsirf to generate Z scores"
        " that are used to determine enriched peptides."
)

# action set up for diffEnrich tsv pipeline
plugin.pipelines.register_function(
    function=diffEnrich_tsv,
    inputs={},
    outputs=diffEnrich_shared_outputs,
    parameters={
        "raw_data_filepath": Str,
        **diffEnrich_tsv_params,
        **diffEnrich_shared_parameters
    },
    input_descriptions=None,
    output_descriptions=None,
    parameter_descriptions={
        "raw_data_filepath": "Raw data matrix in .tsv format.",
        **diffEnrich_tsv_param_descriptions,
        **diffEnrich_shared_parameter_description
    },
    name="diffEnrich tsv Pepsirf Pipeline",
    description="Uses the diff normaization from pepsirf to generate Z scores"
        " that are used to determine enriched peptides."
)

plugin.pipelines.register_function(
    function=diffEnrich_deconv,
    inputs={
        "raw_data": FeatureTable[RawCounts],
        "negative_control": FeatureTable[Normed],
        "bins": PeptideBins,
        "thresh_file": EnrichThresh,
        "linked":Link,
        "id_name_map":PepsirfDMP,
    },
    outputs=[
        ("dir_out", DeconvBatch),
        ("score_per_round", ScorePerRound),
        ("map_dir", PeptideAssignmentMap),
        ("col_sum", FeatureTable[Normed]),
        ("diff", FeatureTable[NormedDifference]),
        ("diff_ratio", FeatureTable[NormedDiffRatio]),
        ("zscore", FeatureTable[Zscore]),
        ("zscore_nan", ZscoreNan),
        ("sample_names", InfoSNPN),
        ("read_counts", InfoSumOfProbes),
        ("rc_boxplot", Visualization),
        ("enrich", PairwiseEnrichment),
        ("enrich_count_boxplot", Visualization),
        ("zscore_scatter", Visualization),
        ("colsum_scatter", Visualization),
        ("zenrich_scatter", Visualization)
    ],
    parameters={
        "deconv_threshold": Int,
        "mapfile_suffix": Str,
        "outfile_suffix": Str,
        "scoring_strategy": Str,
        "score_filtering": Bool,
        "score_tie_threshold": Float,
        "score_overlap_threshold": Float,
        "single_threaded": Bool,
        "remove_file_types": Bool,
        **diffEnrich_shared_parameters,
    },
    input_descriptions=None,
    output_descriptions=None,
    parameter_descriptions={
        **diffEnrich_shared_parameter_description
    },
    name="diffEnrich deconv Pepsirf Pipeline",
    description="Uses the diff normalization from pepsirf to generate z scores"
        " that are used to determine enriched peptides and"
        " **ADD DECONV DESCRIPTION**"
)

plugin.pipelines.register_function(
    function=diffEnrich_deconv_tsv,
    inputs={},
    outputs=[
        ("dir_out", DeconvBatch),
        ("score_per_round", ScorePerRound),
        ("map_dir", PeptideAssignmentMap),
        *diffEnrich_shared_outputs
    ],
    parameters={
        "deconv_threshold": Int,
        "mapfile_suffix": Str,
        "outfile_suffix": Str,
        "scoring_strategy": Str,
        "score_filtering": Bool,
        "score_tie_threshold": Float,
        "score_overlap_threshold": Float,
        "single_threaded": Bool,
        "remove_file_types": Bool,
        "raw_data_tsv": Str,
        "negative_control_tsv": Str,
        "bins_tsv": Str,
        "thresh_file_tsv": Str,
        "linked_tsv": Str,
        "id_name_map_tsv": Str,
        **diffEnrich_shared_parameters,
    },
    input_descriptions=None,
    output_descriptions=None,
    parameter_descriptions={
        **diffEnrich_shared_parameter_description
    },
    name="diffEnrich deconv Pepsirf Pipeline",
    description="Uses the diff normalization from pepsirf to generate z scores"
        " that are used to determine enriched peptides and"
        " **ADD DECONV DESCRIPTION**"
)

# shared outputs for redoDemux and redoDemux tsv pipeline
shared_outputs = [
    ("filtered_counts_output", FeatureTable[RawCounts]),
    ("bad_output_zscores", Visualization), 
    ("good_output_zscores", Visualization)
]

# shared paremters for redoDemux and redoDemux tsv pipeline
shared_parameters = {
    "count_thresh": Float,
    "max_zeros": Int,
    "drop_samp_out": Str,
    "pairs_file": Str,
    "log_normalization": Bool,
    "correlation_threshold": Float,
    "filtered_good_tsv_dir": Str
}

# shared parameter descriptions for redoDemux and redoDemux tsv pipeline
shared_parameter_description = {
    "count_thresh": "Minimum sequence count to not be filtered out. If None is provided,"
        " default is 2x the total number of unique peptides.",
    "max_zeros": "Maximum number of zero counts a sequence needs to not be filtered out. If"
        " None is provided, default is 25% of the total number of unique peptides",
    "drop_samp_out": "Filepath to output sample names that are filtered out of matrix.",
    "pairs_file": "Optional tab-delimited file containing pairs of"
            " sample names. Include if samples are longitudinal.",
    "log_normalization": "Run a log normalization on each of the sets of"
            " scores before running a correlation test on them.",
    "correlation_threshold": "Set a threshold value; anything below the"
        " value will be considered a bad correlation score, and anything"
        " above will be considered a good correlation score.",
    "filtered_good_tsv_dir": "Directory where filtered good matrix files are output to."
}

# shared paremters for redoDemux and redoDemux tsv pipeline
redoDemux_shared_parameters = {
    "seq": Str,
    "read_per_loop": Int,
    "num_threads": Int,
    "phred_base": Int,
    "phred_min_score": Int,
    "sindex": Str,
    "translate_aggregates": Bool,
    "concatemer": Bool,
    "sname": Str,
    "index1": Str,
    "index2": Str,
}

# shared parameter descriptions for redoDemux and redoDemux tsv pipeline
redoDemux_shared_parameter_description = {
    "seq": "Positional information for the DNA tags. This argument must be"
        " passed in the same format specified for 'index1'.",
    "read_per_loop": "The number of fastq records read a time. A higher"
        " value will result in more memory usage by the program, but will"
        " also result in fewer disk accesses, increasing performance of"
        " the program.",
    "num_threads": "Number of threads to use for analyses.",
    "phred_base": "Phred base to use when parsing fastq quality scores."
        " Valid options include 33 or 64.",
    "phred_min_score": "The minimum average phred-scaled quality score for"
        " the DNA tag portion of a read for it to be considered for"
        " matching. This means that if the average phred33/64 score for a"
        " read at the expected locations of the DNA tag is not at least"
        " this then the read will be discarded.",
    "sindex": "Used to specify the header for the index 1 and optional"
        " index 2 column in the samplelist. This is an alternative to"
        " using the '--fif'' option.",
    "translate_aggregates": "Include this flag to use translation-based"
        " aggregation. In this mode, counts for nt sequences will be"
        " combined if they translate into the same aa sequence. Note: When"
        " this mode is used, the name of the aggregate sequence will be"
        " the sequence that was a result of the translation. Therefore,"
        " this mode is most appropriate for use with reference-independent"
        " demultiplexing.",
    "concatemer": "Concatenated adapter/primer sequences (optional). The"
        " presence of this sequence within a read indicates that the"
        " expected DNA tag is not present. If supplied, the number of"
        " times this concatemer is recorded in the input file is"
        " reported.",
    "sname": "Used to specify the header for the sample name column in the"
        " samplelist. By default 'SampleName' is set as the column header"
        " name.",
    "index1": "Positional information for index1 (i.e barcode 1). This"
        " argument must be passed as 3 comma-separated values. The first"
        " item represents the (0-based) expected start position of the"
        " first index; the second represents the length of the first"
        " index; and the third represents the number of mismatches that"
        " are tolerated for this index. An example is '--index1 12,12,1'."
        " This says that the index starts at (0-based) position 12, the"
        " index is 12 nucleotides long, and if a perfect match is not"
        " found, then up to one mismatch will be tolerated.",
    "index2": "Positional information for index2, optional. This argument"
        " must be passed in the same format specified for '--index1'. If"
        " '--input2' is provided, this positional information is assummed"
        " to refer to the reads contained in this second, index-only fastq"
        " file. If '--input_r2' is NOT provided, this positional"
        " information is assumed to refer to the reads contained in the"
        " '--input_r1' fastq file.",
}

# action set up for redoDemux module
plugin.pipelines.register_function(
    function=redoDemux,
    inputs={
        "input_r1": DemuxFastq,
        "input_r2": DemuxFastq,
        "index": DemuxIndex,
        "samplelist": DemuxSampleList,
        "fif": DemuxFif,
        "library": DemuxLibrary,
        **diffEnrich_inputs
    },
    outputs=[
        ("raw_counts_output", FeatureTable[RawCounts]),
        ("diagnostic_output", DemuxDiagnostic),
        *shared_outputs,
        *diffEnrich_good_outputs,
        *diffEnrich_shared_outputs
    ],
    parameters={
        **redoDemux_shared_parameters,
        **shared_parameters,
        **diffEnrich_shared_parameters
    },
    input_descriptions={
       "input_r1": "Fastq-formatted file containing reads with DNA tags. If"
            " PepSIRF was NOT compiled with Zlib support, this file must be"
            " uncompressed. If PepSIRF was compiled with Zlib support, then"
            " this file can be uncompressed or compressed using gzip. In this"
            " case, the file format will be automatically determined.",
        "input_r2": "Optional index-only fastq file. If PepSIRF was NOT"
            " compiled with Zlib support, this file must be uncompressed. If"
            " PepSIRF was compiled with Zlib support, then this file can be"
            " uncompressed or compressed using gzip. In this case, the file"
            " format will be automatically determined. Note that if this"
            " argument is not supplied, only 'index1' will be used to identify"
            " samples.",
        "index": "Name of fasta-formatted file containing forward and"
            " (potentially) reverse index sequences. Sequence names must match"
            " exactly with those supplied in the 'samplelist'.",
        "samplelist": "A tab-delimited list of samples with a header row and"
            " one sample per line. This file must contain at least one index"
            " column and one sample name column. Multiple index columns may be"
            " included. This file can also include additional columns that"
            " will not be used for the demultiplexing. Specify which columns"
            " to use with the '--sname', '--sindex1', and '--sindex2' flags."
            " If '-fif' is used, then only '-sname' will be used.",
        "fif": "The flexible index file can be provided as an alternative to"
            " the '--index1' and '--index2' options. The file must use the"
            " following format: a tab-delimited file with 5 ordered columns:"
            " 1) index name, which should correspond to a header name in the"
            " sample sheet, 2) read name, which should be either 'r1' or 'r2'"
            " (not case-sensitive) to specify whether the index is in"
            " '--input_r1' or '--input_r2', 3) index start location (0-based,"
            " inclusive), 4) index length and 5) number of mismatched to"
            " allow. '--index1', '--index2', '--sname', '--sindex1', and"
            " 'sindex2' will be ignored if this option is provided.",
        "library": "Fasta-formatted file containing reference DNA tags. If"
            " this flag is not included, reference-independent demultiplexing"
            " will be performed. In reference-independent mode, each sequence"
            " in the region specified by '--seq' will be considered its own"
            " reference, and the observed sequences will be used as the row"
            " names in the output count matrix.",
        **diffEnrich_input_descriptions
    },
    output_descriptions=None,
    parameter_descriptions={
        **redoDemux_shared_parameter_description,
        **shared_parameter_description,
        **diffEnrich_shared_parameter_description
    },
    name="redoDemux Pepsirf Pipeline",
    description=""
)

# action set up for redoDemux tsv pipeline
plugin.pipelines.register_function(
    function=redoDemux_tsv,
    inputs={},
    outputs=[
        ("raw_counts_output", FeatureTable[RawCounts]),
        ("diagnostic_output", DemuxDiagnostic),
        *shared_outputs,
        *diffEnrich_good_outputs,
        *diffEnrich_shared_outputs
    ],
    parameters={
        "input_r1_filepath": Str,
        "input_r2_filepath": Str,
        "index_filepath": Str,
        "samplelist_filepath": Str,
        "fif_filepath": Str,
        "library_filepath": Str,
        **redoDemux_shared_parameters,
        **shared_parameters,
        **diffEnrich_shared_parameters,
        **diffEnrich_tsv_params
    },
    input_descriptions=None,
    output_descriptions=None,
    parameter_descriptions={
        "input_r1_filepath": "Filepath to file containing reads with DNA tags. If"
            " PepSIRF was NOT compiled with Zlib support, this file must be"
            " uncompressed. If PepSIRF was compiled with Zlib support, then"
            " this file can be uncompressed or compressed using gzip. In this"
            " case, the file format will be automatically determined.",
        "input_r2_filepath": "Filepath to optional index-only fastq file. If PepSIRF was NOT"
            " compiled with Zlib support, this file must be uncompressed. If"
            " PepSIRF was compiled with Zlib support, then this file can be"
            " uncompressed or compressed using gzip. In this case, the file"
            " format will be automatically determined. Note that if this"
            " argument is not supplied, only 'index1' will be used to identify"
            " samples.",
        "index_filepath": "Filepath to name of fasta-formatted file containing forward and"
            " (potentially) reverse index sequences. Sequence names must match"
            " exactly with those supplied in the 'samplelist'.",
        "samplelist_filepath": "Filepath to a tab-delimited list of samples with a header row and"
            " one sample per line. This file must contain at least one index"
            " column and one sample name column. Multiple index columns may be"
            " included. This file can also include additional columns that"
            " will not be used for the demultiplexing. Specify which columns"
            " to use with the '--sname', '--sindex1', and '--sindex2' flags."
            " If '-fif' is used, then only '-sname' will be used.",
        "fif_filepath": "Filepath to the flexible index file can be provided as an alternative to"
            " the '--index1' and '--index2' options. The file must use the"
            " following format: a tab-delimited file with 5 ordered columns:"
            " 1) index name, which should correspond to a header name in the"
            " sample sheet, 2) read name, which should be either 'r1' or 'r2'"
            " (not case-sensitive) to specify whether the index is in"
            " '--input_r1' or '--input_r2', 3) index start location (0-based,"
            " inclusive), 4) index length and 5) number of mismatched to"
            " allow. '--index1', '--index2', '--sname', '--sindex1', and"
            " 'sindex2' will be ignored if this option is provided.",
        "library_filepath": "Filepath to file containing reference DNA tags. If"
            " this flag is not included, reference-independent demultiplexing"
            " will be performed. In reference-independent mode, each sequence"
            " in the region specified by '--seq' will be considered its own"
            " reference, and the observed sequences will be used as the row"
            " names in the output count matrix.",
        **redoDemux_shared_parameter_description,
        **shared_parameter_description,
        **diffEnrich_shared_parameter_description,
        **diffEnrich_tsv_param_descriptions
    },
    name="redoDemux tsv Pepsirf Pipeline",
    description=""
)

# action set up for redoDemux module
plugin.pipelines.register_function(
    function=redoNoDemux,
    inputs={
        "input_raw_counts": FeatureTable[RawCounts],
        **diffEnrich_inputs
    },
    outputs=[
        *shared_outputs,
        *diffEnrich_good_outputs,
        *diffEnrich_shared_outputs
    ],
    parameters={
        **shared_parameters,
        **diffEnrich_shared_parameters
    },
    input_descriptions={
        "input_raw_counts": "FeatureTable containing raw PepSIRF counts matrix for filtering.",
        **diffEnrich_input_descriptions
    },
    output_descriptions=None,
    parameter_descriptions={
        **shared_parameter_description,
        **diffEnrich_shared_parameter_description
    },
    name="redoNoDemux Pepsirf Pipeline",
    description=""
)

# action set up for redoDemux tsv pipeline
plugin.pipelines.register_function(
    function=redoNoDemux_tsv,
    inputs={},
    outputs=[
        *shared_outputs,
        *diffEnrich_good_outputs,
        *diffEnrich_shared_outputs
    ],
    parameters={
        "input_raw_counts_filepath": Str,
        **shared_parameters,
        **diffEnrich_shared_parameters,
        **diffEnrich_tsv_params
    },
    input_descriptions=None,
    output_descriptions=None,
    parameter_descriptions={
        "input_raw_counts_filepath": "Filepath to .tsv containing raw PepSIRF counts matrix for filtering.",
        **shared_parameter_description,
        **diffEnrich_shared_parameter_description,
        **diffEnrich_tsv_param_descriptions

    },
    name="redoNoDemux tsv Pepsirf Pipeline",
    description=""
)
