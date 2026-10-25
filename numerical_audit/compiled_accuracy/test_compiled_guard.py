"""Regression for a watchdog stop whose systemd-run supervisor returns zero.

The actual premise was observed in i07-preflight-03: CPU stop, SIGTERM and no
completed audit, despite wrapper exit zero. These tests execute the full CLI
control path with owned-service responses mocked; they start no real services.
"""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('guard_under_test_'+HERE.name,HERE/'guarded.py')
guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)

@pytest.mark.parametrize('reason',['CPU','Disk','Wall'])
def test_budget_stop_is_failure_even_when_supervisor_reports_success(tmp_path,monkeypatch,reason):
    logs=tmp_path/'logs';output=tmp_path/'output';output.mkdir();calls=[]
    if reason=='Disk':(output/'retained').write_bytes(b'x'*(1024**2+1))
    class Process:
        stopped=False
        def poll(self):return 0 if self.stopped else None
        def wait(self,**kwargs):return 0
    process=Process()
    def systemctl(args,**kwargs):
        calls.append(args)
        assert args[:2]==['systemctl','--user'] and args[3].startswith('fermionic-future-probe-')
        if args[2]=='stop':process.stopped=True;return SimpleNamespace(returncode=0,stdout='')
        assert args[2]=='show'
        return SimpleNamespace(returncode=0,stdout='CPUUsageNSec='+('2000000000' if reason=='CPU' else '0'))
    ticks=iter([0,2,2,2]) if reason=='Wall' else iter([0,0,0,0])
    monkeypatch.setattr(guard.time,'monotonic',lambda:next(ticks,2 if reason=='Wall' else 0))
    monkeypatch.setattr(guard,'available',lambda:6*1024**3)
    monkeypatch.setattr(guard.subprocess,'Popen',lambda *args,**kwargs:process)
    monkeypatch.setattr(guard.subprocess,'run',systemctl)
    monkeypatch.setattr('sys.argv',['guarded.py','--label','probe','--logs',str(logs),'--cpu-seconds','1',
        '--wall-seconds','1','--disk-mib','1','--output',str(output),'--','python','never-executed.py'])
    with pytest.raises(SystemExit) as exit:guard.main()
    record=json.loads((logs/'probe.resources.json').read_text())
    assert exit.value.code!=0 and record['exit_code']!=0
    assert record['supervisor_exit_code']==0 and record['termination']==reason+' budget exhausted'
    assert any(args[2]=='stop' for args in calls)

@pytest.mark.parametrize('supervisor_status',[0,7])
def test_normal_completion_preserves_actual_supervisor_status(tmp_path,monkeypatch,supervisor_status):
    class Process:
        def poll(self):return supervisor_status
        def wait(self,**kwargs):return supervisor_status
    monkeypatch.setattr(guard,'available',lambda:6*1024**3)
    monkeypatch.setattr(guard.subprocess,'Popen',lambda *args,**kwargs:Process())
    monkeypatch.setattr('sys.argv',['guarded.py','--label','normal','--logs',str(tmp_path),'--cpu-seconds','1',
        '--wall-seconds','1','--disk-mib','1','--output',str(tmp_path/'output'),'--','python','never-executed.py'])
    with pytest.raises(SystemExit) as exit:guard.main()
    record=json.loads((tmp_path/'normal.resources.json').read_text())
    assert exit.value.code==record['exit_code']==record['supervisor_exit_code']==supervisor_status
    assert record['termination'] is None

@pytest.mark.parametrize('reason',['CPU','Disk','Wall'])
def test_completed_job_final_limits_are_enforced(tmp_path,monkeypatch,reason):
    """A child may finish between polls; its retained final usage still counts."""
    logs=tmp_path/'logs';output=tmp_path/'output';output.mkdir()
    class Process:
        def poll(self):return 0
        def wait(self,**kwargs):return 0
    def start(*args,**kwargs):
        if reason=='CPU':
            kwargs['stdout'].write('User time (seconds): 0.80\nSystem time (seconds): 0.30\n')
            kwargs['stdout'].flush()
        if reason=='Disk':(output/'retained').write_bytes(b'x'*(1024**2+1))
        return Process()
    ticks=iter([0,2,2,2]) if reason=='Wall' else iter([0,0,0,0])
    monkeypatch.setattr(guard.time,'monotonic',lambda:next(ticks,2 if reason=='Wall' else 0))
    monkeypatch.setattr(guard,'available',lambda:6*1024**3)
    monkeypatch.setattr(guard.subprocess,'Popen',start)
    monkeypatch.setattr('sys.argv',['guarded.py','--label','final','--logs',str(logs),'--cpu-seconds','1',
        '--wall-seconds','1','--disk-mib','1','--output',str(output),'--','python','never-executed.py'])
    with pytest.raises(SystemExit) as exit:guard.main()
    record=json.loads((logs/'final.resources.json').read_text())
    assert exit.value.code!=0 and record['exit_code']!=0
    assert record['supervisor_exit_code']==0 and record['termination']==reason+' budget exhausted'
