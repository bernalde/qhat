"""Targeted metadata contracts; synthetic rows/costs are not study results.

No Julia child, native objective, chemistry or physical evaluation is invoked.
The coordinator runs this file under its existing finite test resource guard.
"""
import copy
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

HERE=Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module;spec.loader.exec_module(module);return module
r=load('resume_contracts',HERE/'calibration_resume.py')
m=load('resume_prior_cost_contracts',HERE/'calibration_completion.py')
c=load('resume_unchanged_fit_contracts',HERE/'calibration.py')


class ResumeContracts(unittest.TestCase):
    def resources(self):
        return dict(schema='calibration-resume-resources-v1',declared_before_use=True,
            cumulative_pool=dict(CPU_seconds=840,wall_seconds=6000,disk_MiB=256),
            attempt=dict(CPU_seconds=360,wall_seconds=1800,disk_MiB=128),memory_profile=dict(r.MEMORY))
    def costs(self):
        return [dict(charged_CPU_seconds=238.072,wall_seconds=952.422,output_bytes=2*1024**2),
                dict(charged_CPU_seconds=201.857,wall_seconds=2991.473,output_bytes=2*1024**2)]
    def rows(self,indices,commit):
        return [dict(index=i,verified=True,producer=dict(commit=commit),artifact=dict(path='row'+str(i),sha256='a'*64))
                for i in indices]
    def test_missingness_is_exact_and_independent_of_quality(self):
        rows=self.rows(range(13),'023');before=copy.deepcopy(rows)
        reused,missing=r.partition_rows(list(reversed(rows)))
        self.assertEqual([row['index'] for row in reused],list(range(13)))
        self.assertEqual(missing,list(range(13,33)))
        self.assertEqual(rows,before)
        scattered=self.rows((0,3,30),'023')
        self.assertEqual(r.partition_rows(scattered)[1],[i for i in range(33) if i not in (0,3,30)])
        with self.assertRaises(r.Refusal):r.partition_rows(rows+[dict(rows[0],pretend_better_score=0.)])
        with self.assertRaises(r.Refusal):r.partition_rows([dict(rows[0],verified=False)])
        with self.assertRaises(r.Refusal):r.partition_rows([dict(rows[0],index=True)])
    def test_combination_retains_prior_and_resume_producers_without_overwrite(self):
        old=self.rows(range(13),'023');new=self.rows(range(13,33),'resume');before=copy.deepcopy(old)
        combined=r.combine_rows(old,new)
        self.assertEqual(len(combined),33)
        self.assertTrue(all(row['producer']['commit']=='023' for row in combined[:13]))
        self.assertTrue(all(row['producer']['commit']=='resume' for row in combined[13:]))
        self.assertEqual(old,before)
        with self.assertRaises(r.Refusal):r.combine_rows(old,new[:-1])
        with self.assertRaises(r.Refusal):r.combine_rows(old,new+[dict(old[0],producer=dict(commit='resume'))])
    def test_all_interrupted_costs_reduce_remaining_pool(self):
        costs=self.costs();pool,attempt=r.resource_plan(costs,self.resources())
        self.assertAlmostEqual(sum(row['charged_CPU_seconds'] for row in costs),439.929)
        self.assertAlmostEqual(pool['CPU_seconds']-sum(row['charged_CPU_seconds'] for row in costs),400.071)
        self.assertEqual(attempt['CPU_seconds'],360)
        unknown=self.costs();unknown[1]['charged_CPU_seconds']=None
        with self.assertRaises(r.Refusal):r.resource_plan(unknown,self.resources())
        exhausted=self.costs();exhausted[1]['wall_seconds']=4000
        with self.assertRaises(r.Refusal):r.resource_plan(exhausted,self.resources())
        invalid=self.resources();invalid['attempt']['disk_MiB']=math.inf
        with self.assertRaises(r.Refusal):r.resource_plan(costs,invalid)
        old_profile=self.resources();old_profile['memory_profile']['MemoryMax']='512M'
        with self.assertRaises(r.Refusal):r.resource_plan(costs,old_profile)
    def test_unchanged_complete_64_probe_fit_retains_losses_and_floor_refusals(self):
        probes=[]
        for case,parent in (('CO_r1p13',1.),('CH4_r1p09',2.)):
            for i in range(1,33):
                delta=(.1,.2,-.1,0.)[(i-1)%4]
                probes.append(dict(case=case,index=i,parent_J=parent,native_J=parent*(1+delta),verified=True))
        before=copy.deepcopy(probes);fit=c.fit_temperature(probes)
        self.assertEqual(probes,before)
        self.assertEqual(fit['positive_delta_count'],32)
        self.assertAlmostEqual(fit['pooled_positive_median'],.15)
        self.assertAlmostEqual(math.exp(-fit['pooled_positive_median']/fit['theta_start']),.8)
        self.assertAlmostEqual(math.exp(-fit['pooled_positive_median']/fit['theta_end']),.01)
        floor=copy.deepcopy(probes);floor[5]['verified']=False
        with self.assertRaises(c.Refusal):c.fit_temperature(floor)
        with self.assertRaises(c.Refusal):c.fit_temperature(probes[:-1])
    def test_new_receipt_requires_prospective_verification_ram_and_terminal_costs(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary);log=folder/'job.log';log.write_text('CPU time consumed: 2min 1.234s\n')
            path=folder/'receipt.json';command=['python','resume.py','work']
            profile=dict(m.PROFILE,**r.MEMORY)
            receipt=dict(command=command,completed_utc='fixture',settings=profile,
                environment=dict(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1'),
                guard_error=None,cleanup_errors=[],exit_code=0,termination=None,elapsed=500.,log=str(log),
                process_measurements=dict(user_cpu_seconds=120.2,system_cpu_seconds=1.),
                final_accounting=dict(process_CPU_timing_available=True,CPU_seconds=121.2,output_bytes=1024))
            path.write_text(json.dumps(receipt))
            _,cost=r.new_receipt(m,path,command,profile)
            self.assertAlmostEqual(cost['charged_CPU_seconds'],121.235)
            self.assertEqual(cost['wall_seconds'],500.)
            self.assertFalse(cost['failed_attempt'])
            bad=copy.deepcopy(receipt);bad['settings']['MemoryMax']='512M';path.write_text(json.dumps(bad))
            with self.assertRaises(r.Refusal):r.new_receipt(m,path,command,profile)
            unknown=copy.deepcopy(receipt);unknown['final_accounting']['process_CPU_timing_available']=False
            log.write_text('no terminal CPU footer\n');path.write_text(json.dumps(unknown))
            with self.assertRaises(r.Refusal):r.new_receipt(m,path,command,profile)


if __name__=='__main__':unittest.main()
