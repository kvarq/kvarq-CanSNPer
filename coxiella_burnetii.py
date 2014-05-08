
from _CanSNPer import make_CanSNPer_testsuite, VERSION, GENES_COMPATIBILITY

coxiella_burnetii = make_CanSNPer_testsuite(
        VERSION,
        {'AE016828.2': 'AE016828.2.fa'},
        'coxiella_burnetii_snp.txt',
        'coxiella_burnetii_tree.txt'
    )

