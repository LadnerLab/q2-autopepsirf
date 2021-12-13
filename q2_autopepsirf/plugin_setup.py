import importlib

from qiime2.plugin import (Plugin, TypeMap, Str, List,
                            MetadataColumn, Categorical,
                            Int, Range, Visualization)
import q2_autopepsirf

from q2_types.feature_table import FeatureTable
from q2_pepsirf.format_types import (RawCounts, Normed, NormedDifference,
                NormedDiffRatio, PeptideBins, Zscore, ZscoreNan, InfoSNPN,
                EnrichThresh, PairwiseEnrichment, InfoSumOfProbes)
from q2_autopepsirf.actions.diffEnrich import diffEnrich

plugin = Plugin("autopepsirf", version=q2_autopepsirf.__version__,
                website="https://github.com/LadnerLab/q2-autopepsirf")

plugin.pipelines.register_function(
    function=diffEnrich,
    inputs={
        'raw_data': FeatureTable[RawCounts | Normed],
        'negative_control': FeatureTable[RawCounts | Normed],
        'bins': PeptideBins,
        'thresh_file': EnrichThresh
    },
    outputs=[
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
        'negative_id': Str,
        'negative_names': List[Str],
        'pepsirf_binary': Str,
        'exact_z_thresh': Str,
        'raw_constraint': Int % Range(0, None),
        'exact_zenrich_thresh': List[Str],
        'step_z_thresh': Int % Range(1, None),
        'upper_z_thresh': Int % Range(2, None),
        'lower_z_thresh': Int % Range(1, None),
        'pepsirf_tsv_dir': Str,
        'tsv_base_str': Str
    },
    input_descriptions={
        'raw_data': "Raw data matrix.",
        'negative_control': "Name of FeatureTable matrix file containing data for sb samples.",
        'bins': "Name of the file containing bins, one bin per line, as output by the bin module. Each bin contains a "
                "tab-delimited list of peptide names.",
        'thresh_file': "The name of a tab-delimited file containing one tab-delimited matrix filename "
                    "and threshold(s), one per line. If providing more than z score matrix."
    },
    output_descriptions=None,
    parameter_descriptions={
        'negative_id': "Optional approach for identifying negative controls. Provide a unique string at the start of all "
                    "negative control samples.",
        'negative_names': "Optional approach for identifying negative controls. "
                        "Space-separated list of negative control sample names.",
        'pepsirf_binary': "The binary to call pepsirf on your system.",
        'exact_z_thresh': "Individual Exact z score threshold separated by a comma for creation of threshold file"
                        " to run pepsirf's enrich module (Ex: 6,10 or 5,30)",
        'raw_constraint': "The minimum total raw count across all peptides for a sample to be "
                        "included in the analysis.This provides a way to impose a minimum read "
                        "count for a sample to be evaluated.",
        'exact_zenrich_thresh': "List of exact z score thresholds either individual or combined. "
                            "List MUST BE in descending order. (Example argument: '--p-exact-zenrich-thresh 25 10 3' "
                            "or '--p-exact-zenrich-thresh 6,25 4,10 1,3')",
        "step_z_thresh": "Integar to increment z-score thresholds.",
        "upper_z_thresh": "Upper limit of z-score thresholds (non-inclusive).",
        "lower_z_thresh": "Lower limit of z-score thresholds (inclusive).",
        "pepsirf_tsv_dir": "Provide a directory path. Must also provide tsv-base-str",
        "tsv_base_str": "The base name for the output tsv files excluding ay extensions, typcally the raw data filename "
                        "(EX: --p-tsv-base-str raw_data). Must also provide pepsirf-tsv-dir"
    },
    name='diffEnrich Pepsirf Pipeline',
    description="Uses the diff normaization from "
                "pepsirf to generate Z scores that are used to determine enriched peptides"
)
