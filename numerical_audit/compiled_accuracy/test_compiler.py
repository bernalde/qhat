"""Independent literal-unitary tests for the unexecuted prototype compiler."""
from fractions import Fraction
import itertools
import gzip
import json
import numpy as np
import pytest
from compiler import Compiler,Gate,read_trace,resources,write_trace

I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)
H=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
S=np.diag([1,1j])

def embed(matrix,q,n):
    out=np.array([[1.+0j]])
    for at in reversed(range(n)):out=np.kron(out,matrix if at==q else I)
    return out

def literal_gate(g,n):
    if g.kind=='cx':
        # Independent basis permutation, including nonlocal CX for its target.
        out=np.zeros((1<<n,1<<n),complex)
        for col in range(1<<n):out[col^((1<<g.b) if (col>>g.a)&1 else 0),col]=1
        return out
    if g.kind=='rz':
        a=float(g.angle)/2;matrix=np.diag([np.exp(-1j*a),np.exp(1j*a)])
    else:matrix={'h':H,'s':S,'sdg':S.conj().T}[g.kind]
    return embed(matrix,g.a,n)

def matrix(gates,n):
    out=np.eye(1<<n,dtype=complex)
    for g in gates:out=literal_gate(g,n)@out
    return out

def pauli(word,n):
    out=np.eye(1<<n,dtype=complex)
    for q,p in word:out=embed({'X':X,'Y':Y,'Z':Z}[p],q,n)@out
    return out

WORDS=[[(q,p) for q,p in enumerate(labels) if p!='I'] for labels in itertools.product('IXYZ',repeat=3) if labels!=('I','I','I')]

@pytest.mark.parametrize('word',WORDS)
@pytest.mark.parametrize('angle',[Fraction(2,7),Fraction(-5,11)])
def test_complete_pauli_unitary_including_Y_and_disconnected_support(word,angle):
    c=Compiler(3);c.rotation(word,angle);P=pauli(word,3)
    expected=np.cos(float(angle))*np.eye(8)-1j*np.sin(float(angle))*P
    np.testing.assert_allclose(matrix(c.gates,3),expected,atol=3e-14,rtol=0)
    assert all(g.kind!='cx' or abs(g.a-g.b)==1 for g in c.gates)

@pytest.mark.parametrize('control,target',[(a,b) for a in range(4) for b in range(4) if a!=b])
def test_restored_routing_is_full_nonlocal_CX(control,target):
    c=Compiler(4);c.remote_cx(control,target)
    np.testing.assert_array_equal(matrix(c.gates,4),literal_gate(Gate('cx',control,target),4))
    assert len(c.gates)==6*(abs(control-target)-1)+1

@pytest.mark.parametrize('electrons',range(5))
def test_HF_preparation_has_exact_phase(electrons):
    c=Compiler(4);c.prepare_hf(electrons);a=matrix(c.gates,4)[:,0]
    expected=np.zeros(16,complex);expected[(1<<electrons)-1]=1
    np.testing.assert_allclose(a,expected,atol=2e-15,rtol=0)
    assert resources(4,c.gates)['two_qubit_gates']==0

def test_exact_adjacent_rewrites_preserve_full_unitary():
    raw=[Gate('h',0),Gate('h',0),Gate('s',1),Gate('sdg',1),Gate('cx',0,1),Gate('cx',0,1),
         Gate('rz',0,angle=Fraction(1,7)),Gate('rz',0,angle=Fraction(2,7)),Gate('rz',0,angle=Fraction(-3,7))]
    c=Compiler(2)
    for g in raw:c.add(g)
    assert c.gates==[]
    np.testing.assert_allclose(matrix(raw,2),np.eye(4),atol=2e-15,rtol=0)
    # No cancellation/commuting across an intervening disjoint gate.
    for g in [Gate('h',0),Gate('h',1),Gate('h',0)]:c.add(g)
    assert len(c.gates)==3

def test_depth_counts_include_single_qubit_dependency_and_parallel_layers():
    gates=[Gate('h',0),Gate('s',2),Gate('cx',0,1),Gate('rz',1,angle=Fraction(1,3)),Gate('cx',1,2)]
    assert resources(3,gates)==dict(gates=5,single_qubit_gates=3,two_qubit_gates=2,unit_duration_ASAP_depth=4,CX_dependency_depth=2)

def test_complete_exact_trace_roundtrip_and_deterministic_archive(tmp_path):
    c=Compiler(3);c.prepare_hf(2);c.rotation([(0,'Y'),(2,'X')],Fraction(-13,37))
    a,b=tmp_path/'a.gz',tmp_path/'b.gz';identity={'source':'fixture','steps':3}
    assert write_trace(a,3,c.gates,identity)==write_trace(b,3,c.gates,identity)
    iterator=read_trace(a);header=next(iterator);assert header['identity']==identity
    restored=list(iterator);assert restored==c.gates and resources(3,restored)==resources(3,c.gates)
    np.testing.assert_array_equal(matrix(restored,3),matrix(c.gates,3))

def test_reject_invalid_register_support_and_primitive_route():
    for n in (0,9,True):
        with pytest.raises(ValueError):Compiler(n)
    c=Compiler(3)
    for word in ([],[(0,'X'),(0,'Z')],[(3,'X')],[(0,'A')]):
        with pytest.raises(ValueError):c.rotation(word,Fraction(1,7))
    with pytest.raises(ValueError):c.add(Gate('cx',0,2))
    with pytest.raises(ValueError):c.rotation([(0,'X')],.1)

def test_reject_corrupt_trace_qubit_and_zero_rotation(tmp_path):
    for payload in (['cx',0,2],['rz',0,'0','1'],['h',True]):
        p=tmp_path/(str(len(list(tmp_path.iterdir())))+'.gz')
        with gzip.open(p,'wt') as f:
            f.write(json.dumps({'compiler':__import__('compiler').VERSION,'n':3})+'\n'+json.dumps(payload)+'\n')
        with pytest.raises(ValueError):list(read_trace(p))
