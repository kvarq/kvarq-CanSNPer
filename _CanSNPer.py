
VERSION = '1.0'
GENES_COMPATIBILITY = '0.1'

from kvarq.genes import Genome, Genotype, Test, Reference, SNP, Testsuite
from kvarq.genes import TemplateFromGenome, Gene
from kvarq.log import lo

import os.path
from functools import reduce

class SNPNode:
    '''
    tree node with for SNP based hierarchical genotyping
    '''

    def __init__(self, parent, name, tests):
        self.parent = parent
        self.name = name
        self.tests = tests
        self.children = list()

        if parent:
            self.parent.children.append(self)

    def get(self, child_name):
        for child in self.children:
            if child.name == child_name:
                return child
        return None

    def path(self):
        if self.parent is None:
            return self.name
        return self.parent.path() + ';' + self.name

    def is_ancestor(self, candidate):
        if self.parent:
            if self.parent == candidate:
                return True
            return self.parent.is_ancestor(candidate)
        return False

class SNPTreeTestsuite(Testsuite):

    def __init__(self, tests, root, version):
        super(SNPTreeTestsuite, self).__init__(tests, version)
        self.root = root

    def walk_tree(self, ret, node, coverages):
        '''
        yields a list of tree-nodes (leaves first) that match the coverages
        '''
        for child in node.children:
            n = 0
            for test in child.tests:
                coverage = coverages[test]
                if test.template.validate(coverage):
                    n += 1
            if n == len(child.tests):
                self.walk_tree(ret, child, coverages)
                ret.append(child)
            elif n>0:
                lo.debug('discarding SNPNode "%s" : %d<%d' % (
                        child.path(), n, len(child.tests)))

    def _analyse(self, coverages):
        nodes = list()
        self.walk_tree(nodes, self.root, coverages)
        i = 0
        # discard all but longest of common paths
        while i < len(nodes):
            j = i + 1
            while j < len(nodes):
                if nodes[i].is_ancestor(nodes[j]):
                    lo.debug('pruning "%s" < "%s"' % (nodes[j].path(), nodes[i].path()))
                    del nodes[j]
                else:
                    j +=1
            i += 1
        return [node.path() for node in nodes]


def make_CanSNPer_testsuite(version, genome_fnames, snp_fname, tree_fname, **kwargs):
    '''
    generate a testsuite that classifies an organism based on a
    hierarchically defined list of single nucleotide polymorphisms
    '''
    base = os.path.dirname(__file__)
    genomes = {}
    for genome, fname in genome_fnames.items():
        genomes[genome] = Genome(os.path.join(base, fname))

    # create SNPs
    tests = dict()
    refs = dict()
    organism = None
    strain = None
    for line in file(os.path.join(base, snp_fname)).readlines():
        if line[0] == '#': continue

        parts = line.rstrip('\n\r').split('\t')

        name = parts[0]
        assert name not in tests, 'SNP "%s" defined twice in "%s"' % (parts[0], snp_fname)

        if organism is None:
            organism = parts[1]
        assert organism == parts[1], "all SNPs must be of same organism"
        ref = refs.setdefault(parts[2], Reference(parts[2]))

        strain = parts[3]
        pos = int(parts[4])
        newbase = parts[5]
        oldbase = parts[6]

        if newbase == '/':
            newbases = set('ACGT') - set(oldbase)
        else:
            newbases = set(newbase)

        tests[name] = list()
        for newbase in newbases:

            snp = SNP(genome=genomes[strain], pos=pos, base=newbase, orig=oldbase, force=True)
            genotype = Genotype(name)
            tests[name].append(Test(template=snp, genotype=genotype, reference=ref))

    # construct Test-tree
    root = None
    for line in file(os.path.join(base, tree_fname)).readlines():
        if root is None:
            assert ';' not in line, 'first line in file "%s" must be root' % tree_fname
            root = SNPNode(None, line.rstrip('\n\r'), [])
            continue
        node = root
        for part in line.rstrip('\n\r').split(';')[1:]:
            assert part in tests, 'SNP "%s" not found in "%s"' % (part, snp_fname)
            x = node.get(part)
            if x is None:
                x = SNPNode(node, name=part, tests=tests[part])
            node = x

    testlist = reduce(lambda x,y: x+y, tests.values())
    return SNPTreeTestsuite(testlist, root, version, **kwargs)

