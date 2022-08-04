from q2_autopepsirf.actions.diffEnrich import diffEnrich
from q2_autopepsirf.actions.diffEnrich_tsv import diffEnrich_tsv
from q2_autopepsirf.actions.diffEnrich_deconv import diffEnrich_deconv
from q2_autopepsirf.actions.diffEnrich_deconv_tsv import diffEnrich_deconv_tsv

__all__ = ['diffEnrich', 'diffEnrich_tsv', 'diffEnrich_deconv', 'diffEnrich_deconv_tsv']
from . import _version
__version__ = _version.get_versions()['version']
