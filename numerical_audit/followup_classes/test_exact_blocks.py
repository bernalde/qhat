import itertools
import math
import unittest
import numpy as np
from openfermion import FermionOperator, QubitOperator, get_sparse_operator, hermitian_conjugated
from exact_blocks import (exact_parents, decompose, pauli_commutator_term, commutator,
    internally_commuting, number_commutator, graph_and_color, suffix_greedy, suffix_objective)


def matrix(block,n=4):
    op=QubitOperator(); op.terms={k:complex(v) for k,v in block.items()}
    return get_sparse_operator(op,n_qubits=n).toarray()


class ExactBlockTests(unittest.TestCase):
    def test_pauli_products_exhaustive_two_modes(self):
        words=[tuple((q,p) for q,p in enumerate(labels) if p!='I')
               for labels in itertools.product('IXYZ',repeat=2)]
        mats={w:matrix({w:1},2) for w in words}
        for a,b in itertools.product(words,repeat=2):
            dest,scale=pauli_commutator_term(a,b)
            actual=mats[a]@mats[b]-mats[b]@mats[a]
            expected=np.zeros((4,4),complex) if scale==0 else 1j*scale*mats[dest]
            np.testing.assert_array_equal(actual,expected)

    def test_class_commutation_uses_sum_cancellation(self):
        number={((0,'Z'),):-1,((1,'Z'),):-1}
        hop={((0,'X'),(1,'X')):1,((0,'Y'),(1,'Y')):1}
        self.assertTrue(any(pauli_commutator_term(a,b)[1] for a in number for b in hop))
        self.assertEqual(commutator(number,hop),{})
        _,graph,_=graph_and_color([number,hop])
        self.assertEqual(graph['edge_count'],0)
        self.assertEqual(graph['color_count'],1)

    def test_exact_parent_and_ownership_conserve(self):
        raw=FermionOperator('0^ 1',.25)+FermionOperator('1^ 0',.25)
        raw+=FermionOperator('0^ 0 1^ 1',.375)
        parents,total,den,info=exact_parents(raw,0)
        for blocks in decompose(parents,total)[:3]:
            for block in blocks:
                self.assertTrue(internally_commuting(block))
                self.assertFalse(number_commutator(block))
        target=get_sparse_operator(raw,n_qubits=2).toarray()
        represented=matrix({k:v/den for k,v in total.items()},2)
        shift=np.trace(target-represented)/4
        np.testing.assert_array_equal(target,represented+shift*np.eye(4))
        self.assertEqual(info['canonical_hermitian_adjustment_fermionic_l1'],0)

    def test_greedy_objective_and_strict_swaps(self):
        blocks=[{((0,'X'),):3},{((0,'Y'),):2},{((0,'Z'),):1},{((0,'X'),):-2}]
        _,_,cs=graph_and_color(blocks,True)
        order,info=suffix_greedy(len(blocks),cs)
        self.assertEqual(sorted(order),list(range(4)))
        self.assertLessEqual(int(info['final_objective_integer']),int(info['greedy_objective_integer']))
        for order in itertools.permutations(range(4)):
            manual=0
            for i in range(4):
                suffix={}
                for j in order[i+1:]:
                    for k,v in blocks[j].items():suffix[k]=suffix.get(k,0)+v
                manual+=sum(abs(v) for v in commutator(blocks[order[i]],suffix).values())
            self.assertEqual(suffix_objective(order,cs),manual)

    def test_parent_averaging_is_explicit(self):
        raw=FermionOperator('0^ 1',.25)+FermionOperator('1^ 0',.5)
        parents,_,den,info=exact_parents(raw,0)
        self.assertEqual(abs(parents[0]['weight']/den),.375)
        self.assertEqual(info['canonical_hermitian_adjustment_fermionic_l1'],.25)

    def test_four_orbital_classes_and_physical_ratio(self):
        normalized={}
        for perm in itertools.permutations(range(4)):
            a,b,c,d=perm
            m=FermionOperator(((a,1),(b,1),(c,0),(d,0)),1)
            raw=m+hermitian_conjugated(m)
            parents,_,den,_=exact_parents(raw,0)
            self.assertEqual(len(parents),1)
            p=parents[0]
            normalized.setdefault(p['components'],[]).append(perm)
            self.assertEqual(len(p['coefficients']),8)
        self.assertEqual(sorted(map(len,normalized.values())),[8,8,8])
        for weights in ((1,3,2),(0,2,2),(1,1,0)):
            raw=FermionOperator()
            for perm,w in zip(((0,1,2,3),(0,2,1,3),(0,3,1,2)),weights):
                a,b,c,d=perm;m=FermionOperator(((a,1),(b,1),(c,0),(d,0)),w)
                raw+=m+hermitian_conjugated(m)
            parents,total,_,_=exact_parents(raw,0)
            self.assertEqual(sum(len(p['coefficients']) for p in parents),4*len(total))

if __name__=='__main__': unittest.main()
