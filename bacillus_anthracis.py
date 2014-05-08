
from _CanSNPer import make_CanSNPer_testsuite, VERSION, GENES_COMPATIBILITY

bacillus_anthracis = make_CanSNPer_testsuite(
        VERSION,
        dict(Ames='Ames.fa'),
        'bacillus_anthracis_snp.txt',
        'bacillus_anthracis_tree.txt'
    )

