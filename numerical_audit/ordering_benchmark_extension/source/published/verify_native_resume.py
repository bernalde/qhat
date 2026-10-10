"""Resume an actual molecular multiprecision V2 search from a partial checkpoint."""
import hashlib,json,resource,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main():
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));base=HERE/'v2/native/f2_coeff1p00';payload=(base/'input.txt').read_text();gold=(base/'accepted.tsv').read_text().splitlines();half=len(gold)//2
    work=HERE/'resume_check';work.mkdir(exist_ok=True);checkpoint=work/'accepted.tsv';checkpoint.write_text('\n'.join(gold[:half])+'\n');started=time.monotonic()
    result=subprocess.run([str(HERE/'v2_exact'),str(checkpoint)],input=payload,text=True,capture_output=True,check=True)
    recovered=checkpoint.read_text().splitlines()
    if recovered!=gold:raise RuntimeError('Partial native resume differs from uninterrupted exact history')
    d=json.loads((HERE/'v2/designs/f2_coeff1p00.json').read_text());expected=d['methods']['tranter_v2_l1_first_order']['metadata']['eq4_order']
    if json.loads(result.stdout)!=expected:raise RuntimeError('Partial native resume final order mismatch')
    report={'case':'f2_coeff1p00','total_terms':len(gold),'accepted_terms_replayed':half,'remaining_terms_recomputed':len(gold)-half,'all_exact_history_norms_identical':True,'final_order_identical':True,'seconds':time.monotonic()-started,'input_sha256':hashlib.sha256(payload.encode()).hexdigest(),'checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((HERE/'v2_exact.cpp').read_bytes()).hexdigest()}
    (HERE/'native_partial_resume_check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
