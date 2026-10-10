#!/usr/bin/env python3
"""Resume-friendly completion supervisor; never alters frozen plan selection."""
import json,os,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent));from common import CASES

def main():
    start=time.monotonic()
    while not all((HERE/'v2/designs'/f'{case}.json').exists() for case in CASES):
        if time.monotonic()-start>6*3600:raise RuntimeError('Native V2 preparation did not complete within6hours; inspect checkpoint/log before resuming')
        time.sleep(15)
    for command in [['run.py','freeze','v2'],['check_scaling.py'],['run.py','run','v2'],['run.py','collect','v2'],['summarize.py']]:
        subprocess.run([sys.executable,str(HERE/command[0]),*command[1:]],check=True)
    print('V2 CAMPAIGN COMPLETE',flush=True)
if __name__=='__main__':main()
