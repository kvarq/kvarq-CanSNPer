
from _CanSNPer import make_CanSNPer_testsuite, VERSION, GENES_COMPATIBILITY

francisella_tularensis = make_CanSNPer_testsuite(
        VERSION,
        {
            'OSU18': 'OSU18.fa',
            'SCHUS4.1': 'SCHUS4.1.fa',
            'SCHUS4.2': 'SCHUS4.2.fa'
        },
        'francisella_tularensis_snp.txt',
        'francisella_tularensis_tree.txt'
    )

