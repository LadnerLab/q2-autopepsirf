import importlib

from qiime2.plugin import (Plugin, TypeMap, Str, List,
                            MetadataColumn, Categorical,
                            Int, Range, Visualization)

from q2_types.feature_table import FeatureTable
from q2_pepsirf.format_types import (RawCounts, Normed, NormedDifference,
                NormedDiffRatio, PeptideBins, Zscore, ZscoreNan, InfoSNPN,
                EnrichThresh, PairwiseEnrichment, InfoSumOfProbes)
from q2_autopepsirf.actions.diffEnrich import diffEnrich

plugin = Plugin("autopepsirf", version="0.0.1.dev",
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
        ("zscore_out", FeatureTable[Zscore]),
        ("nan_out", ZscoreNan),
        ("sample_names", InfoSNPN),
        ("read_counts", InfoSumOfProbes),
        ("rc_boxplot", Visualization),
        ("enrich_out", PairwiseEnrichment),
        ("enrich_count_boxplot", Visualization),
        ("zscore_scatter", Visualization),
        ("colsum_scatter", Visualization)
    ],
    parameters={
        'negative_id': Str,
        'negative_names': List[Str],
        'pepsirf_binary': Str,
        'exact_z_thresh': Str,
        'raw_constraint': Int % Range(0, None)
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
        'exact_z_thresh': "Exact z score thresholds either individual or combined.",
        'raw_constraint': "The minimum total raw count across all peptides for a sample to be "
                        "included in the analysis.This provides a way to impose a minimum read "
                        "count for a sample to be evaluated.",
    },
    name='diffEnrich Pepsirf Pipeline',
    description="Uses the diff normaization from "
                "pepsirf to generate Z scores that are used to determine enriched peptides"
)
