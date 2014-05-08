
import unittest

from kvarq.log import set_debug
set_debug()

import os.path
import sys
sys.path.insert(0, os.path.dirname(__file__))

from bacillus_anthracis import bacillus_anthracis
from yersinia_pestis import yersinia_pestis
from coxiella_burnetii import coxiella_burnetii 

def create_coverages(tests, snps):
    coverages = dict([(test, False) for test in tests])
    for test in tests:
        if test.genotype.identifier in snps:
            coverages[test] = True
    return coverages

class GenesCanSNPer(unittest.TestCase):

    def test_bacillus_anthracis(self):

        # patch tests
        for test in bacillus_anthracis.tests:
            test.template.validate = lambda x: x

        # child / parent nodes must not be included
        coverages = create_coverages(bacillus_anthracis.tests, 'A/B.Br.001;B.Br.003'.split(';'))
        ret = bacillus_anthracis._analyse(coverages)
        assert 'A/B.Br.001;B.Br.003' in ret
        assert 'A/B.Br.001;B.Br.003;B.Br.004' not in ret
        assert 'A/B.Br.001' not in ret

        # missing node
        coverages = create_coverages(bacillus_anthracis.tests, 'A/B.Br.001;A.Br.006;A.Br.004;A.Br.003;A_Br_002;A.Br.001'.split(';'))
        ret = bacillus_anthracis._analyse(coverages)
        assert ret == ['A/B.Br.001;A.Br.006;A.Br.004;A.Br.003']

    def test_yersinia_pestis(self):
        # patch tests
        for test in yersinia_pestis.tests:
            test.template.validate = lambda x: x
        # test
        coverages = create_coverages(yersinia_pestis.tests, 'Root;0.PE3.a;III;VI;VII;0.ANT3.a;3.ANT.a;VIII;2.MED3.a;IX;X;2.MED1.a;2.MED1.b;2.MED1.d'.split(';'))
        ret = yersinia_pestis._analyse(coverages)
        assert ret == ['Root;0.PE3.a;III;VI;VII;0.ANT3.a;3.ANT.a;VIII;2.MED3.a;IX;X;2.MED1.a;2.MED1.b;2.MED1.d']


    def test_coxiella_burnetii(self):
        # patch tests
        for test in coxiella_burnetii.tests:
            test.template.validate = lambda x: x
        # test
        coverages = create_coverages(coxiella_burnetii.tests, 'Root;C.7;C.9;C.10'.split(';'))
        ret = coxiella_burnetii._analyse(coverages)
        assert ret == ['Root;C.7;C.9;C.10']

if __name__ == '__main__': unittest.main()
