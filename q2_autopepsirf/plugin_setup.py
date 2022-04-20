import importlib

from qiime2.plugin import (Plugin, TypeMap, Str, List,
                            MetadataColumn, Categorical,
                            Int, Range, Visualization, Float,
                            Bool)
import q2_autopepsirf

from q2_types.feature_table import FeatureTable
from q2_pepsirf.format_types import (RawCounts, Normed, NormedDifference,
                NormedDiffRatio, PeptideBins, Zscore, ZscoreNan, InfoSNPN,
                EnrichThresh, PairwiseEnrichment, InfoSumOfProbes, DeconvBatch,
                PeptideAssignmentMap, ScorePerRound, Link, PepsirfDMP)
from q2_autopepsirf.actions.diffEnrich import diffEnrich
from q2_autopepsirf.actions.diffEnrich_tsv import diffEnrich_tsv
from q2_autopepsirf.actions.diffEnrich_deconv import diffEnrich_deconv
from q2_autopepsirf.actions.diffEnrich_deconv_tsv import diffEnrich_deconv_tsv

# This is the plugin object. It is what the framework will load and what an
# interface will interact with. Basically every registration we perform will
# involve this object in some way.
plugin = Plugin("autopepsirf", version=q2_autopepsirf.__version__,
                website="https://github.com/LadnerLab/q2-autopepsirf",
                description="Qiime2 plugin used for the automation of q2-pepsirf and q2-ps-plot.")

# shared outputs for diffEnrich and diffEnrich tsv pipeline
shared_outputs = [
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

# shared paremters for diffEnrich and diffEnrich tsv pipeline
shared_parameters = {
        'negative_id': Str,
        'negative_names': List[Str],
        'pepsirf_binary': Str,
        'exact_z_thresh': Str,
        'exact_cs_thresh': Str,
        'raw_constraint': Int % Range(0, None),
        'exact_zenrich_thresh': List[Str],
        'step_z_thresh': Int % Range(1, None),
        'upper_z_thresh': Int % Range(2, None),
        'lower_z_thresh': Int % Range(1, None),
        'pepsirf_tsv_dir': Str,
        'tsv_base_str': Str,
        'hdi': Float % Range(0.0, 1.0),
        'infer_pairs_source': Bool,
        'flexible_reps_source': Bool,
        's_enrich_source': Bool,
        'user_defined_source': MetadataColumn[Categorical]

    }

# shared parameter descriptions for diffEnrich and diffEnrich tsv pipeline
shared_parameter_description = {
        'negative_id': "Optional approach for identifying negative controls. Provide a unique string at the start of all "
                    "negative control samples.",
        'negative_names': "Optional approach for identifying negative controls. "
                        "Space-separated list of negative control sample names.",
        'pepsirf_binary': "The binary to call pepsirf on your system.",
        'exact_z_thresh': "Individual Exact z score threshold separated by a comma for creation of threshold file"
                        " to run pepsirf's enrich module (Ex: 6,10 or 30)",
        'exact_cs_thresh': "Individual Exact col-sum threshold separated by a comma for creation of threshold file"
                        " to run pepsirf's enrich module (Ex: 6,10 or 30)",
        'raw_constraint': "The minimum total raw count across all peptides for a sample to be "
                        "included in the analysis.This provides a way to impose a minimum read "
                        "count for a sample to be evaluated.",
        'exact_zenrich_thresh': "List of exact z score thresholds either individual or combined. "
                            "List MUST BE in descending order. (Example argument: '--p-exact-zenrich-thresh 25 10 3' "
                            "or '--p-exact-zenrich-thresh 6,25 4,10 1,3')",
        "step_z_thresh": "Integar to increment z-score thresholds.",
        "upper_z_thresh": "Upper limit of z-score thresholds (non-inclusive).",
        "lower_z_thresh": "Lower limit of z-score thresholds (inclusive).",
        "pepsirf_tsv_dir": "Provide a directory path. Must also provide tsv-base-str for output of tsv verison of qza files."
                        " The source_samples file and png boxplot outputs will always be put within this directory.",
        "tsv_base_str": "The base name for the output tsv files excluding ay extensions, typcally the raw data filename "
                        "(EX: --p-tsv-base-str raw_data). Must also provide pepsirf-tsv-dir, if pepsirf-tsv-dir provided without "
                        "tsv-base-str, the default will be 'aps-output'.",
        "hdi": "Alternative approach for discarding outliers prior to calculating mean and stdev. If provided, this "
            "argument will override --trim, which trims evenly from both sides of the distribution. For --hdi, the "
            "user should provide the high density interval to be used for calculation of mean and stdev. For "
            "example, '--hdi 0.95' would instruct the program to utilize the 95% highest density interval (from each "
            "bin) for these calculations.",
        "infer_pairs_source": "Infer sample pairs from names. This option assumes names of replicates will be identical "
                            "with the exception of a final string denoted with a '_'. For example, these names would be "
                            "considered two replicates of the same sample: VW_100_1X_A and VW_100_1X_B",
        "flexible_reps_source": "will infer the number of replicates for each sample based on sample names, and will not "
                            "require any specific number of replicates for inclusion. Therefore, some samples may have a "
                            "single replicate, some may have 2, 3, 4 etc. And all replicates of a given sample will be "
                            "considered for determining enriched peptides.",
        "s_enrich_source": "All samples will be processed individually as samples with only one replicate",
        "user_defined_source": "Metadata file containing all sample names and their source groups. "
                            "Used to create pairs tsv to run pepsirf enrich module."
    }

# action set up for diffEnrich module
plugin.pipelines.register_function(
    function=diffEnrich,
    inputs={
        'raw_data': FeatureTable[RawCounts],
        'negative_control': FeatureTable[Normed],
        'bins': PeptideBins,
        'thresh_file': EnrichThresh
    },
    outputs=shared_outputs,
    parameters=shared_parameters,
    input_descriptions={
        'raw_data': "Raw data matrix.",
        'negative_control': "Name of FeatureTable matrix file containing data for sb samples.",
        'bins': "Name of the file containing bins, one bin per line, as output by the bin module. Each bin contains a "
                "tab-delimited list of peptide names.",
        'thresh_file': "The name of a tab-delimited file containing one tab-delimited matrix filename "
                    "and threshold(s), one per line. If providing more than z score matrix."
    },
    output_descriptions=None,
    parameter_descriptions=shared_parameter_description,
    name='diffEnrich Pepsirf Pipeline',
    description="Uses the diff normaization from "
                "pepsirf to generate Z scores that are used to determine enriched peptides"
)

# action set up for diffEnrich tsv pipeline
plugin.pipelines.register_function(
    function=diffEnrich_tsv,
    inputs={},
    outputs=shared_outputs,
    parameters={
        'raw_data_tsv': Str,
        'negative_control_tsv': Str,
        'bins_tsv': Str,
        'thresh_file_tsv': Str,
        **shared_parameters
    },
    input_descriptions=None,
    output_descriptions=None,
    parameter_descriptions={
        'raw_data_tsv': "Raw data matrix in .tsv format.",
        'negative_control_tsv': "Name of .tsv matrix file containing data for sb samples.",
        'bins_tsv': "Name of the file containing bins, one bin per line, as output by the bin module. Each bin contains a "
                        "tab-delimited list of peptide names.",
        'thresh_file_tsv': "The name of a tab-delimited file containing one tab-delimited matrix filename "
                            "and threshold(s), one per line. If providing more than z score matrix.",
        **shared_parameter_description
    },
    name='diffEnrich tsv Pepsirf Pipeline',
    description="Uses the diff normaization from "
                "pepsirf to generate Z scores that are used to determine enriched peptides"
)

plugin.pipelines.register_function(
    function=diffEnrich_deconv,
    inputs={
        'raw_data': FeatureTable[RawCounts],
        'negative_control': FeatureTable[Normed],
        'bins': PeptideBins,
        'thresh_file': EnrichThresh,
        'linked':Link,
        'id_name_map':PepsirfDMP,
    },
    outputs={
        ('dir_out', DeconvBatch),
        ('score_per_round', ScorePerRound),
        ('map_dir', PeptideAssignmentMap),
        *shared_outputs
    },
    parameters={
        'deconv_threshold': Int,
        'mapfile_suffix' : Str,
        'outfile_suffix' : Str,
        'scoring_strategy' : Str,
        'score_filtering' : Bool,
        'score_tie_threshold' : Float,
        'score_overlap_threshold' : Float,
        'single_threaded' : Bool,
        'remove_file_types' : Bool,
        **shared_parameters,
    },
    input_descriptions=None,
    output_descriptions=None,
    parameter_descriptions={
        **shared_parameter_description
    },
    name='diffEnrich deconv Pepsirf Pipeline',
    description="Uses the diff normalization from "
                "pepsirf to generate z scores that are used to determine enriched peptides"
                "and **ADD DECONV DESCRIPTION**"
)

plugin.pipelines.register_function(
    function=diffEnrich_deconv_tsv,
    inputs={
        
    },
    outputs={
        ('dir_out', DeconvBatch),
        ('score_per_round', ScorePerRound),
        ('map_dir', PeptideAssignmentMap),
        *shared_outputs
    },
    parameters={
        'deconv_threshold': Int,
        'mapfile_suffix' : Str,
        'outfile_suffix' : Str,
        'scoring_strategy' : Str,
        'score_filtering' : Bool,
        'score_tie_threshold' : Float,
        'score_overlap_threshold' : Float,
        'single_threaded' : Bool,
        'remove_file_types' : Bool,
        'raw_data_tsv': Str,
        'negative_control_tsv': Str,
        'bins_tsv': Str,
        'thresh_file_tsv': Str,
        'linked_tsv' : Str,
        'id_name_map_tsv' : Str,
        **shared_parameters,
    },
    input_descriptions=None,
    output_descriptions=None,
    parameter_descriptions={
        **shared_parameter_description
    },
    name='diffEnrich deconv Pepsirf Pipeline',
    description="Uses the diff normalization from "
                "pepsirf to generate z scores that are used to determine enriched peptides"
                "and **ADD DECONV DESCRIPTION**"
)
