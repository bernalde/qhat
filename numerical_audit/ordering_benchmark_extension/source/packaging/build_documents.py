#!/usr/bin/env python3
"""Compile both manuscript entry points with durable logs and explicit paths."""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('--paper',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--entry',choices=['main','theory_review','theory_report'],required=True);a=p.parse_args();root=a.paper.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True);entry={'main':'main.tex','theory_review':'theory/theory_review.tex','theory_report':'theory_report.tex'}[a.entry];steps=[]
 commands=[['pdflatex','-no-shell-escape','-interaction=nonstopmode','-halt-on-error','-recorder','-output-directory='+str(out),entry],['bibtex',a.entry]]
 commands+= [commands[0],commands[0]]
 for i,cmd in enumerate(commands):
  env=os.environ.copy();env['BIBINPUTS']=str(root)+os.pathsep+env.get('BIBINPUTS','');cwd=out if cmd[0]=='bibtex' else root
  start=time.monotonic();r=subprocess.run(cmd,cwd=cwd,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);log=out/f'command_{i}.log';log.write_text(r.stdout);steps.append({'command':cmd,'returncode':r.returncode,'seconds':time.monotonic()-start})
  if r.returncode:raise RuntimeError('Document build failed: '+str(log)+'\n'+r.stdout[-5000:])
 pdf=out/(a.entry+'.pdf');log=(out/(a.entry+'.log')).read_text();warnings=[line for line in log.splitlines() if 'Warning' in line or 'Overfull' in line or 'Underfull' in line];record={'entry':entry,'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'pdf_bytes':pdf.stat().st_size,'steps':steps,'warnings':warnings};(out/'build_record.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'pdf':str(pdf),'warnings':warnings}))
if __name__=='__main__':main()
