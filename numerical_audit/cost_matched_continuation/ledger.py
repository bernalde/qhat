"""Archive the complete scientific cutoff before analysis; preserve all jobs."""
import argparse
import datetime
import common as k
from construct import copy_evidence


def archive(root, campaign):
    k.require(not (root/'campaign/LEDGER_CUTOFF.json').exists(), 'Preserve original ledger cutoff')
    copy_evidence(campaign/'RESOURCE_POOLS.json', root/'campaign/RESOURCE_POOLS.json')
    pools = k.read(campaign/'RESOURCE_POOLS.json')
    guard = campaign/pools['guard']['path']
    k.require(k.sha(guard) == pools['guard']['sha256'], 'Frozen guard bytes differ')
    copy_evidence(guard, root/pools['guard']['path'])
    copy_evidence(campaign/'campaign/run_guarded_job.py', root/'campaign/run_guarded_job.py')
    entries = []
    for plan_path in sorted((campaign/'jobs').glob('*.PLAN.json')):
        plan = k.read(plan_path)
        label = plan['label']
        entry = dict(label=label, phase=plan['phase'])
        for kind, name in (('plan', 'jobs/'+label+'.PLAN.json'), ('final', 'jobs/'+label+'.FINAL.json'),
                           ('receipt', 'logs/'+label+'.resources.json')):
            source = campaign/name
            k.require(source.is_file(), 'Unfinished attempt blocks scientific cutoff: '+label)
            copy_evidence(source, root/'campaign'/name)
            entry[kind+'_path'] = 'campaign/'+name
            entry[kind+'_sha256'] = k.sha(source)
        copy_evidence(campaign/'logs'/(label+'.log'), root/'campaign/logs'/(label+'.log'))
        entries.append(entry)
    k.write_new(root/'campaign/LEDGER_CUTOFF.json', dict(schema='cost-matched-ledger-v1',
        cutoff_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        resource_pools=dict(path='campaign/RESOURCE_POOLS.json', sha256=k.sha(campaign/'RESOURCE_POOLS.json')),
        jobs=entries, scope='All scientific jobs through cutoff, including failures; later analysis/report/review/packaging recorded separately'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=k.Path, required=True)
    parser.add_argument('--campaign', type=k.Path, required=True)
    args = parser.parse_args()
    archive(args.root, args.campaign)
