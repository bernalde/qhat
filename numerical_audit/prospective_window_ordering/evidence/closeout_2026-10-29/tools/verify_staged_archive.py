"""Check explicit archive inventory, staged Git bytes and high-signal secret markers."""
import argparse, hashlib, json, re, subprocess
from pathlib import Path

PREFIX = 'numerical_audit/prospective_window_ordering/'
ARCHIVE = PREFIX + 'evidence/run_2026-10-28/'
MAX_BYTES = 256 * 1024**2
SECRET = re.compile(rb'(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,}|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----)')

def need(ok, message):
    if not ok: raise ValueError(message)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('worktree', type=Path)
    p.add_argument('manifest_sha256')
    p.add_argument('report_sha256')
    p.add_argument('output', type=Path)
    a = p.parse_args()
    manifest_path = a.worktree / ARCHIVE / 'COPY_MANIFEST.json'
    manifest_bytes = manifest_path.read_bytes()
    need(hashlib.sha256(manifest_bytes).hexdigest() == a.manifest_sha256, 'copy manifest changed')
    manifest = json.loads(manifest_bytes)
    need(manifest['status'] == 'complete' and len(manifest['files']) <= 16382, 'copy manifest status/count')
    expected = {ARCHIVE + r['destination']: (r['sha256'], r['bytes']) for r in manifest['files']}
    need(len(expected) == len(manifest['files']), 'duplicate archive destination')
    expected[ARCHIVE + 'COPY_MANIFEST.json'] = (a.manifest_sha256, len(manifest_bytes))
    report = a.worktree / PREFIX / 'REPORT.md'
    expected[PREFIX + 'REPORT.md'] = (a.report_sha256, report.stat().st_size)
    actual = {ARCHIVE + str(v.relative_to(a.worktree / ARCHIVE)) for v in (a.worktree / ARCHIVE).rglob('*') if v.is_file()}
    need(actual == {n for n in expected if n.startswith(ARCHIVE)}, 'filesystem archive inventory differs')
    staged = subprocess.check_output(['git', 'diff', '--cached', '--name-only', '-z'], cwd=a.worktree).split(b'\0')
    need({n.decode() for n in staged if n} == set(expected), 'staged change inventory differs')
    entries = subprocess.check_output(['git', 'ls-files', '--stage', '-z'], cwd=a.worktree)
    objects = {}
    for entry in entries.split(b'\0'):
        if not entry: continue
        metadata, name = entry.split(b'\t', 1)
        mode, oid, stage = metadata.split()
        path = name.decode()
        if path in expected:
            need(stage == b'0' and mode == b'100644', 'unexpected staged mode/stage')
            objects[path] = oid.decode()
    need(set(objects) == set(expected), 'missing staged objects')
    result, total = [], 0
    child = subprocess.Popen(['git', 'cat-file', '--batch'], cwd=a.worktree, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    try:
        for name in sorted(expected):
            sha, size = expected[name]
            total += size
            need(total <= MAX_BYTES, 'bounded staged byte scan exceeded')
            child.stdin.write((objects[name] + '\n').encode()); child.stdin.flush()
            header = child.stdout.readline().split()
            need(len(header) == 3 and header[1] == b'blob' and int(header[2]) == size, 'staged object header/size')
            raw = child.stdout.read(size)
            need(child.stdout.read(1) == b'\n' and len(raw) == size, 'staged object framing')
            need(hashlib.sha256(raw).hexdigest() == sha, 'staged bytes differ: ' + name)
            need(not SECRET.search(raw), 'high-signal credential marker in selected archive')
            result.append(dict(path=name, sha256=sha, bytes=size, git_blob=objects[name]))
        child.stdin.close()
        need(child.wait(timeout=10) == 0, 'Git object reader failed')
    finally:
        if child.poll() is None: child.kill(); child.wait()
    links = []
    for target in re.findall(r'\]\(([^)]+)\)', report.read_text()):
        if ':' not in target and not target.startswith('#'):
            need((report.parent / target.split('#', 1)[0]).is_file(), 'missing report link: ' + target)
            links.append(target)
    output = dict(schema='staged-archive-git-byte-check-v1', status='passed', files=len(result), bytes=total,
                  manifest_sha256=a.manifest_sha256, report_sha256=a.report_sha256,
                  report_file_links_checked=links, high_signal_secret_markers=0, inventory=result,
                  limitation='Byte/publication validation only; disclosed historical source gaps remain.')
    with a.output.open('x') as f: json.dump(output, f, indent=2, allow_nan=False); f.write('\n')
    print('passed', len(result), 'staged files;', total, 'bytes;', len(links), 'local report links')

if __name__ == '__main__':
    main()
