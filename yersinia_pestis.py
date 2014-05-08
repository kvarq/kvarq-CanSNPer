
from _CanSNPer import make_CanSNPer_testsuite, VERSION, GENES_COMPATIBILITY

yersinia_pestis = make_CanSNPer_testsuite(
        VERSION,
        dict(CO92='CO92.fa'),
        'yersinia_pestis_snp.txt',
        'yersinia_pestis_tree.txt'
    )

