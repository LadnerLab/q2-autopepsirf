from q2_autopepsirf.actions.diffEnrich import diffEnrich
from q2_autopepsirf.actions.diffEnrich_tsv import diffEnrich_tsv

__all__ = ['diffEnrich', 'diffEnrich_tsv']
from . import _version
__version__ = _version.get_versions()['version']
