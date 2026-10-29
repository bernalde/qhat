"""Independent metadata checker for the appended postreview CPU sensitivity view."""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import stat

ORIGINAL_SHA = "db2e4271d53fc140168352cf70a430d4cfecfe8b939db0ed2d7738477cf387f0"
SCHEMA = "cpu-pareto-cost-sensitivity-v1"
STAGE = "postreview_numerical_cost_sensitivity_not_primary"
DOMINANCE = "costA<=costB+1e-12 and errorA<=errorB+1e-12; costA<costB-1e-12 or errorA<errorB-1e-12. Postreview CPU-cost sensitivity only."
TOL, QTOL = 1e-12, Fraction(1, 10**12)
MAX_FILE, MAX_TOKENS, MAX_DEPTH = 32 * 1024**2, 1_000_000, 32
CASES = ("LiH_r1p487_631g_14m4e", "LiH_r2p113_631g_14m4e",
         "BeH2_r1p369_631g_14m4e", "BeH2_r1p907_631g_14m4e")
EXPECTED_KEYS = set(itertools.product(CASES, (1, 2), (.4, 1.6), (8, 24)))

def need(condition, reason):
    if not condition:
        raise ValueError(reason)

def safe_path(value, fresh=False):
    supplied = Path(value)
    need(".." not in supplied.parts and "\0" not in str(value), "unsafe path")
    path = supplied.absolute()
    for parent in (path, *path.parents):
        need(not parent.is_symlink(), "symlink file or ancestor")
    if fresh:
        need(path.parent.is_dir() and not path.exists(), "fresh output with existing parent required")
    return path

def read_bounded(value):
    path = safe_path(value)
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, "rb") as stream:
        before = os.fstat(stream.fileno())
        need(stat.S_ISREG(before.st_mode) and before.st_size <= MAX_FILE, "regular input / 32MiB input limit")
        raw = stream.read(MAX_FILE + 1)
        after = os.fstat(stream.fileno())
    need(len(raw) == before.st_size <= MAX_FILE and
         (before.st_size, before.st_mtime_ns, before.st_ctime_ns) ==
         (after.st_size, after.st_mtime_ns, after.st_ctime_ns), "input changed or exceeded limit")
    return raw

def bounded_structure(raw):
    stack, quoted, escaped, tokens = [], False, False, 0
    for byte in raw:
        if quoted:
            if escaped: escaped = False
            elif byte == 92: escaped = True
            elif byte == 34: quoted = False
        elif byte == 34:
            quoted = True
        elif byte in (123, 91):
            stack.append(byte); tokens += 1
            need(len(stack) <= MAX_DEPTH, "JSON depth limit")
        elif byte in (125, 93):
            need(stack and stack.pop() == (123 if byte == 125 else 91), "JSON bracket types")
        elif byte in (44, 58):
            tokens += 1
        need(tokens <= MAX_TOKENS, "JSON structural-token limit")
    need(not stack and not quoted and not escaped, "unfinished JSON structure")

def unique_object(items):
    result = {}
    for key, value in items:
        need(key not in result, "duplicate JSON key")
        result[key] = value
    return result

def finite_float(text):
    value = float(text)
    need(math.isfinite(value), "nonfinite JSON number")
    return value

def forbidden_constant(text):
    raise ValueError("nonfinite JSON constant")

def decode(raw):
    bounded_structure(raw)
    return json.loads(raw, object_pairs_hook=unique_object, parse_float=finite_float,
                      parse_constant=forbidden_constant)

def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)

def unchanged(actual, expected, reason):
    need(canonical(actual) == canonical(expected), reason)

def numeric(value):
    need(type(value) in (int, float), "numeric scalar, excluding bool")
    result = float(value)
    need(math.isfinite(result) and result >= 0, "finite nonnegative cost/error")
    if type(value) is int:
        need(int(result) == value, "integer scalar not exactly representable as binary64")
    return result

def point_numbers(points):
    need(type(points) is list and len(points) == 11, "all eleven original CPU points required")
    result, paths = [], set()
    for point in points:
        need(type(point) is dict and set(point) ==
             {"path", "cost", "projective_error", "physical_floor", "fully_validated", "objective_status"}, "point schema")
        name = point["path"]
        need(isinstance(name, str) and 0 < len(name) <= 256 and
             all(ord(c) >= 32 for c in name) and name not in paths, "unique bounded point path")
        paths.add(name)
        need(type(point["physical_floor"]) is bool and type(point["fully_validated"]) is bool, "point flags")
        status = point["objective_status"]
        need(status is None or (isinstance(status, str) and len(status) <= 128), "objective status type")
        cost, error = numeric(point["cost"]), numeric(point["projective_error"])
        result.append((name, cost, error, Fraction.from_float(cost), Fraction.from_float(error)))
    return result

def float_dominates(a, b):
    return a[1] <= b[1] + TOL and a[2] <= b[2] + TOL and (a[1] < b[1] - TOL or a[2] < b[2] - TOL)

def rational_dominates(a, b):
    return a[3] <= b[3] + QTOL and a[4] <= b[4] + QTOL and (a[3] < b[3] - QTOL or a[4] < b[4] - QTOL)

def relations(points, predicate):
    domination = {b[0]: [a[0] for a in points if a[0] != b[0] and predicate(a, b)] for b in points}
    frontier = [point[0] for point in points if not domination[point[0]]]
    return domination, frontier

def audit(original, sensitivity, report):
    need(type(original) is dict and original.get("mode") == "main", "original main summary")
    groups = original.get("pareto_memberships")
    need(type(groups) is list and len(groups) == 32, "original thirty-two groups")
    need(type(sensitivity) is dict and set(sensitivity) ==
         {"schema", "source_summary_sha256", "stage", "cost_tolerance_seconds", "error_tolerance",
          "groups", "changed_CPU_frontier_indices", "producer_sha256", "scope"}, "sensitivity top-level schema")
    need(sensitivity["schema"] == SCHEMA and sensitivity["stage"] == STAGE and
         sensitivity["source_summary_sha256"] == ORIGINAL_SHA, "sensitivity source/stage binding")
    need(isinstance(sensitivity["producer_sha256"], str) and len(sensitivity["producer_sha256"]) == 64 and
         all(c in "0123456789abcdef" for c in sensitivity["producer_sha256"]), "producer source digest spelling")
    need(isinstance(sensitivity["scope"], str) and 0 < len(sensitivity["scope"]) <= 4096, "bounded scope header")
    report.update(declared_producer_sha256=sensitivity["producer_sha256"], producer_source_binding_checked=False)
    for field in ("cost_tolerance_seconds", "error_tolerance"):
        need(type(sensitivity[field]) is float and sensitivity[field].hex() == TOL.hex(), "declared binary64 tolerance")
    supplied = sensitivity["groups"]
    need(type(supplied) is list and len(supplied) == 32, "sensitivity thirty-two groups")
    seen, changed, exact_changed, differences = set(), [], [], []
    for index, (primary, view) in enumerate(zip(groups, supplied)):
        need(type(primary) is dict and primary.get("expected_paths") == 11, "original expected paths")
        key = tuple(primary[field] for field in ("case", "repetition", "time", "budget_multiplier"))
        need(type(primary["repetition"]) is int and type(primary["budget_multiplier"]) is int and
             key in EXPECTED_KEYS and key not in seen, "original group identity")
        seen.add(key)
        need(type(view) is dict and set(view) ==
             {"index", "case", "repetition", "time", "budget_multiplier", "CPU_seconds",
              "expected_paths", "original_CPU_pareto_paths", "wall_seconds_primary_unchanged", "CPU_frontier_changed"}, "view group schema")
        need(type(view["index"]) is int and view["index"] == index, "zero-based original group index")
        need(type(view["expected_paths"]) is int and view["expected_paths"] == primary["expected_paths"], "expected-path count changed")
        for field in ("case", "repetition", "time", "budget_multiplier"):
            unchanged(view[field], primary[field], "group identity changed")
        cpu, supplied_cpu = primary["sets"]["CPU_seconds"], view["CPU_seconds"]
        need(type(supplied_cpu) is dict and set(supplied_cpu) ==
             {"points", "pareto_paths", "dominated_by", "excluded", "dominance"}, "CPU set schema")
        unchanged(supplied_cpu["points"], cpu["points"], "CPU point data/order/status changed")
        unchanged(supplied_cpu["excluded"], cpu["excluded"], "CPU exclusions changed")
        need(cpu["excluded"] == [] and supplied_cpu["dominance"] == DOMINANCE, "all-point eligibility / predicate")
        unchanged(view["original_CPU_pareto_paths"], cpu["pareto_paths"], "primary CPU frontier changed")
        unchanged(view["wall_seconds_primary_unchanged"], primary["sets"]["wall_seconds"], "primary wall result changed")
        points = point_numbers(cpu["points"])
        floating, float_front = relations(points, float_dominates)
        exact, exact_front = relations(points, rational_dominates)
        unchanged(supplied_cpu["dominated_by"], floating, "producer CPU domination differs from float predicate")
        unchanged(supplied_cpu["pareto_paths"], float_front, "producer CPU frontier differs from float predicate")
        differs = set(float_front) != set(cpu["pareto_paths"])
        need(type(view["CPU_frontier_changed"]) is bool and view["CPU_frontier_changed"] is differs, "changed-group flag")
        if differs: changed.append(index)
        if set(exact_front) != set(cpu["pareto_paths"]): exact_changed.append(index)
        for b in points:
            for a in points:
                if a[0] != b[0] and ((a[0] in floating[b[0]]) != (a[0] in exact[b[0]])):
                    differences.append(dict(index=index, dominator=a[0], dominated=b[0],
                        float_predicate=a[0] in floating[b[0]], exact_rational_predicate=a[0] in exact[b[0]]))
        report["groups_checked"] += 1
        report["points_checked"] += len(points)
    need(seen == EXPECTED_KEYS, "complete original group-key inventory")
    inventory = sensitivity["changed_CPU_frontier_indices"]
    need(type(inventory) is list and all(type(i) is int for i in inventory), "changed inventory integer list")
    unchanged(inventory, changed, "changed-group inventory incomplete or selected")
    report.update(complete=True, producer_float_view_verified=True,
                  exact_rational_agreement=not differences, threshold_rounding_difference=bool(differences),
                  threshold_rounding_differences=differences, changed_CPU_frontier_indices=changed,
                  exact_rational_changed_CPU_frontier_indices=exact_changed, wall_groups_unchanged=32,
                  ordered_nonself_pairs_checked=32 * 11 * 10, objective_status_filter_applied=False)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original_summary")
    parser.add_argument("sensitivity_output")
    parser.add_argument("fresh_check_output")
    args = parser.parse_args()
    output = safe_path(args.fresh_check_output, fresh=True)
    report = dict(schema="cpu-pareto-cost-sensitivity-check-v1", status="incomplete", complete=False,
                  verified=False, groups_checked=0, points_checked=0, original_summary_sha256=ORIGINAL_SHA,
                  rational_tolerances="cost=error=1/1000000000000; input scalars exact Fraction.from_float",
                  scope="Postreview CPU-cost sensitivity only; original primary study remains immutable.")
    try:
        raw = read_bounded(args.original_summary)
        need(hashlib.sha256(raw).hexdigest() == ORIGINAL_SHA, "original summary SHA-256 differs")
        original = decode(raw)
        supplied = read_bounded(args.sensitivity_output)
        report["sensitivity_sha256"] = hashlib.sha256(supplied).hexdigest()
        audit(original, decode(supplied), report)
        report["verified"] = report["exact_rational_agreement"]
        report["status"] = "verified" if report["verified"] else "threshold_rounding_difference"
    except (OSError, ValueError, TypeError, KeyError, IndexError, RecursionError, OverflowError) as error:
        report.update(status="incomplete_or_invalid", failure_type=type(error).__name__, failure=str(error)[:500])
    body = (json.dumps(report, indent=2, allow_nan=False) + "\n").encode()
    need(len(body) <= MAX_FILE, "check output size limit")
    with output.open("xb") as stream: stream.write(body)
    print(report["status"], "groups", report["groups_checked"], "points", report["points_checked"])
    return 0 if report["verified"] else 2

if __name__ == "__main__":
    raise SystemExit(main())
