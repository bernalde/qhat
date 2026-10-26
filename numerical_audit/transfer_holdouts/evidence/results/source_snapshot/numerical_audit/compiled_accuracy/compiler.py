"""Exact-angle prototype Pauli compiler for a fixed nearest-neighbor line.

No synthesis optimality, native-device calibration or hardware execution claim.
The version and source digest must be frozen before any resource experiment.
"""
from fractions import Fraction
import gzip,hashlib,json,time
from pathlib import Path
from typing import NamedTuple

VERSION='line-pauli-exact-angle-0.1.0'
MAX_REGISTER=8
MAX_STORED_GATES=2_000_000
MAX_GENERATED_GATES=12_000_000
ZERO=Fraction(0)

def require(ok,message):
    if not ok:raise ValueError(message)

class Gate(NamedTuple):
    kind:str
    a:int
    b:int=-1
    angle:Fraction=ZERO

class CompilationLimit(ValueError):
    """A declared per-candidate gate limit, distinct from invalid input."""

class Compiler:
    """Only globally adjacent inverses/Rz combination; no commuting rewrites."""
    def __init__(self,n):
        require(type(n) is int and 1<=n<=MAX_REGISTER,'Register bound before allocation')
        self.n=n;self.gates=[];self.generated=0

    def add(self,gate):
        require(gate.kind in ('h','s','sdg','rz','cx'),'Unknown primitive gate')
        require(type(gate.a) is int and 0<=gate.a<self.n,'Gate outside register')
        if gate.kind=='cx':
            require(type(gate.b) is int and 0<=gate.b<self.n and abs(gate.a-gate.b)==1,'Primitive CX must be on the line')
        else:require(gate.b==-1,'Unexpected second qubit')
        require(isinstance(gate.angle,Fraction),'Exact rational angles required')
        if gate.kind!='rz':require(gate.angle==0,'Only Rz has a variable angle')
        if self.generated>=MAX_GENERATED_GATES:raise CompilationLimit('Generated-gate ceiling')
        self.generated+=1
        if gate.kind=='rz' and not gate.angle:return
        if self.gates:
            old=self.gates[-1]
            inverse=((old.kind,gate.kind) in (('h','h'),('s','sdg'),('sdg','s'),('cx','cx'))
                and (old.a,old.b)==(gate.a,gate.b))
            if inverse:self.gates.pop();return
            if old.kind==gate.kind=='rz' and old.a==gate.a:
                self.gates.pop();angle=old.angle+gate.angle
                if angle:self.gates.append(Gate('rz',gate.a,angle=angle))
                return
        if len(self.gates)>=MAX_STORED_GATES:raise CompilationLimit('Stored-gate ceiling')
        self.gates.append(gate)

    def swap(self,a,b):
        require(abs(a-b)==1,'Adjacent SWAP required')
        for c,t in ((a,b),(b,a),(a,b)):self.add(Gate('cx',c,t))

    def remote_cx(self,control,target):
        require(type(control) is int and type(target) is int and 0<=control<self.n and 0<=target<self.n and control!=target,'Invalid routed CX')
        direction=1 if target>control else -1;moved=[];at=control
        while abs(at-target)>1:
            pair=(at,at+direction);self.swap(*pair);moved.append(pair);at+=direction
        self.add(Gate('cx',at,target))
        for pair in reversed(moved):self.swap(*pair)

    def prepare_hf(self,electrons):
        require(type(electrons) is int and 0<=electrons<=self.n,'Invalid HF electron count')
        # H S S H = X exactly, including global phase. No extra gate type.
        for q in range(electrons):
            for kind in ('h','s','s','h'):self.add(Gate(kind,q))

    def rotation(self,word,theta):
        require(isinstance(theta,Fraction),'Exact supplied angle required')
        require(bool(word) and len(word)<=self.n,'Nonempty scalar-removed Pauli word required')
        require(all(type(q) is int and 0<=q<self.n and p in ('X','Y','Z') for q,p in word),'Invalid Pauli support')
        require(len({q for q,p in word})==len(word),'Repeated Pauli support')
        ordered=sorted(word);support=[q for q,p in ordered]
        for q,p in ordered:
            if p=='Y':self.add(Gate('sdg',q))
            if p in ('X','Y'):self.add(Gate('h',q))
        for a,b in zip(support,support[1:]):self.remote_cx(a,b)
        self.add(Gate('rz',support[-1],angle=2*theta))
        for a,b in reversed(list(zip(support,support[1:]))):self.remote_cx(a,b)
        for q,p in reversed(ordered):
            if p in ('X','Y'):self.add(Gate('h',q))
            if p=='Y':self.add(Gate('s',q))

def resources(n,gates):
    require(type(n) is int and 1<=n<=MAX_REGISTER,'Resource register bound')
    layers=[0]*n;cx_layers=[0]*n;count=one=two=0
    for g in gates:
        count+=1;require(count<=MAX_STORED_GATES,'Resource stream bound')
        if g.kind=='cx':
            require(abs(g.a-g.b)==1 and 0<=g.a<n and 0<=g.b<n,'Resource CX outside line')
            d=max(layers[g.a],layers[g.b])+1;layers[g.a]=layers[g.b]=d
            c=max(cx_layers[g.a],cx_layers[g.b])+1;cx_layers[g.a]=cx_layers[g.b]=c;two+=1
        else:
            require(g.kind in ('h','s','sdg','rz') and 0<=g.a<n and g.b==-1,'Resource single-qubit gate invalid')
            layers[g.a]+=1;one+=1
    return dict(gates=count,single_qubit_gates=one,two_qubit_gates=two,
        unit_duration_ASAP_depth=max(layers),CX_dependency_depth=max(cx_layers))

def encode(g):
    if g.kind=='cx':return ['cx',g.a,g.b]
    if g.kind=='rz':return ['rz',g.a,str(g.angle.numerator),str(g.angle.denominator)]
    return [g.kind,g.a]

def write_trace(path,n,gates,identity,*,measurements=None):
    """Deterministic gzip NDJSON retains every exact primitive/angle in order."""
    began=time.process_time();wall=time.monotonic()
    with Path(path).open('xb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as zipped:
            header=dict(compiler=VERSION,n=n,identity=identity,scope='Ideal prototype compiler; no hardware measurement')
            zipped.write((json.dumps(header,sort_keys=True,separators=(',',':'))+'\n').encode())
            for g in gates:zipped.write((json.dumps(encode(g),separators=(',',':'))+'\n').encode())
    serialization=dict(CPU_seconds=time.process_time()-began,wall_seconds=time.monotonic()-wall)
    began=time.process_time();wall=time.monotonic();h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    if measurements is not None:measurements.update(serialization=serialization,hashing=dict(CPU_seconds=time.process_time()-began,wall_seconds=time.monotonic()-wall))
    return h.hexdigest()

def read_trace(path):
    """Stream header and validated primitive tuples without loading a full trace."""
    with gzip.open(path,'rt') as stream:
        header=json.loads(next(stream));require(header['compiler']==VERSION,'Different compiler version')
        n=header['n'];require(type(n) is int and 1<=n<=MAX_REGISTER,'Trace register bound')
        yield header
        for index,line in enumerate(stream):
            require(index<MAX_STORED_GATES,'Trace gate bound');r=json.loads(line)
            require(isinstance(r,list) and r and r[0] in ('cx','h','s','sdg','rz'),'Invalid trace record')
            if r[0]=='cx':
                require(len(r)==3 and all(type(q) is int and 0<=q<n for q in r[1:]) and abs(r[1]-r[2])==1,'Invalid trace CX')
                yield Gate('cx',r[1],r[2])
            elif r[0]=='rz':
                require(len(r)==4 and type(r[1]) is int and 0<=r[1]<n and isinstance(r[2],str) and isinstance(r[3],str) and int(r[3])>0,'Invalid trace Rz')
                angle=Fraction(int(r[2]),int(r[3]));require(angle!=0,'Unreduced zero Rz')
                yield Gate('rz',r[1],angle=angle)
            else:
                require(len(r)==2 and type(r[1]) is int and 0<=r[1]<n,'Invalid trace Clifford')
                yield Gate(r[0],r[1])
