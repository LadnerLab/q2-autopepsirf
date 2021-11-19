import importlib

from qiime2.plugin import Plugin

from q2_types.feature_table import FeatureTable
from q2_pepsirf.format_types import RawCounts
from q2_autopepsirf.actions.diffEnrich import diffEnrich

plugin = Plugin("autopepsirf", version="0.0.1.dev",
                website="https://github.com/LadnerLab/q2-autopepsirf")

plugin.pipelines.register_function(
    function=diffEnrich,
    inputs={
        'raw_data': FeatureTable[RawCounts]
    },
    outputs=[
        ("raw_out", FeatureTable[RawCounts])
    ],
    parameters=None,
    input_descriptions={
        'raw_data': ""
    },
    output_descriptions=None,
    parameter_descriptions=None,
    name='diffEnrich Pepsirf Pipeline',
    description=""
)
