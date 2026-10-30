"""Execution admission for this fixed experiment, separate from valid incumbents.

The parent intentionally has a 3,000-evaluation ceiling. The stronger full-order
control must score four starts and stop normally or at its CPU allowance;
algorithmic evaluation, support, or numerical stops qualify that comparison.
Window work can retain a checked parent at declared budget/capacity stops, but a
scored parent does not turn a runtime, reply, or numerical failure into success.

This pure metadata predicate does not establish guard, cost, hash, objective, or
physical validity. Those remain separate gates. In particular, Julia child exit
zero and an improving selected order are not universal admission requirements.
"""
from __future__ import annotations

import math

PARENT_EVALUATIONS = 3000
CHILD_STOPS = frozenset(('child_cpu_limit', 'child_wall_limit', 'child_evaluation_limit',
                         'child_action_limit', 'child_atom_limit'))
CHILD_REFUSALS = frozenset(value.removeprefix('child_') for value in CHILD_STOPS)
PARENT_BUDGET_STOPS = frozenset(('whole_call_limit', 'evaluation_limit',
    'insufficient_evaluation_budget', 'verification_reserve_unavailable'))
CAPACITY_REASONS = frozenset((
    'reachable vector exceeds max_entries', 'cached reference columns exceed max_entries',
    'reachable unit column exceeds max_entries', 'combined cubic columns exceed max_entries',
    'updated cubic columns exceed max_entries', 'window rational exceeds 4096-bit admission cap',
    'window request exceeds 8 MiB'))
RUNTIME_FAILURES = frozenset(('julia_unavailable', 'unsupported_toml_reader',
    'unsupported_owned_cpu_accounting_platform', 'unsupported_owned_cpu_reaper',
    'owned_child_cpu_unavailable'))
DOMAIN_REFUSALS = frozenset(('unsupported_register_or_unit_cap', 'unsupported_word_cap',
    'unsupported_coherent_reference', 'unsupported_storage_cap', 'window_unavailable',
    'unsupported_empty_or_zero_word_unit', 'unsupported_complex_or_mixed_flip_unit'))
REPLY_FAILURES = frozenset(('invalid_worker_reply', 'missing_worker_reply', 'native_julia_score_mismatch'))


def _integer(value, minimum=0):
    return type(value) is int and value >= minimum


def _order(value):
    return (isinstance(value, (list, tuple)) and bool(value) and
            all(type(item) is int for item in value) and sorted(value) == list(range(len(value))))


def _finite(value):
    if type(value) not in (int, float):
        return False
    try:
        return math.isfinite(value) and value >= 0
    except OverflowError:
        return False


def _scored(result):
    if not isinstance(result, dict) or not _order(result.get('order')):
        return False
    score = result.get('score')
    return (isinstance(score, dict) and _integer(score.get('steps'), 1) and
            all(_finite(score.get(key)) for key in ('value', 'cubic_response')))


def _checkpoint(stats):
    """Recognize API metadata set only after its bounded reply parser succeeds."""
    counters = stats.get('child_counters')
    completed = stats.get('child_completed_evaluations')
    allocated = stats.get('child_allocated_evaluations')
    state, updates = stats.get('child_checkpoint_status'), stats.get('child_checked_updates')
    return (stats.get('backend_runtime_assets_verified') is True and
        type(stats.get('native_parent_revalidation_evaluations')) is int and
        stats['native_parent_revalidation_evaluations'] == 1 and
        stats.get('child_evaluation_scope') == 'checkpoint_prefix' and
        _integer(completed, 1) and _integer(allocated, 1) and completed <= allocated and
        isinstance(counters, dict) and all(_integer(counters.get(key)) for key in
            ('evaluation', 'full_score', 'flat_candidate', 'exact_witness', 'action', 'atom', 'insertion_move')) and
        counters['evaluation'] == completed and counters['full_score'] >= 1 and
        sum(counters[key] for key in ('full_score', 'flat_candidate', 'exact_witness')) <= completed and
        state in ('parent_checked', 'checked_update', 'complete') and _integer(updates) and
        counters['full_score'] >= updates+1 and
        (state != 'parent_checked' or updates == 0) and (state != 'checked_update' or updates >= 1))


def execution_admission(arm, result, parent_order=None, refusal=None):
    """Return a JSON-safe execution classification without changing the result."""
    stats = result.get('stats', {}) if isinstance(result, dict) else {}
    stats = stats if isinstance(stats, dict) else {}
    ws, co, cr = (stats.get(key) for key in ('window_status', 'child_outcome', 'child_refusal'))

    def outcome(admitted, category, reason):
        return dict(version=1, admitted=admitted, category=category, reason=reason,
            window_status=ws if isinstance(ws, str) else None,
            child_outcome=co if isinstance(co, str) else None,
            child_refusal=cr if isinstance(cr, str) else None)

    def deny(category, reason):
        return outcome(False, category, reason)

    def allow(category, reason):
        return outcome(True, category, reason)

    if arm not in ('parent', 'full', 'window'):
        return deny('unknown_or_inconsistent', 'unknown_arm')
    if refusal is not None:
        return deny('worker_refusal', 'outer_worker_refusal')
    status = result.get('status') if isinstance(result, dict) else None
    if status == 'numerical_limit':
        return deny('numerical_refusal', 'numerical_limit')
    if not _scored(result):
        return deny('unscored', 'missing_or_invalid_scored_incumbent')
    if arm == 'parent':
        if status in ('complete', 'cpu_limit'):
            return allow('completed' if status == 'complete' else 'budget_stop', status)
        # The native incremental path reserves its final full-score evaluation.
        if status == 'evaluation_limit' and type(result.get('evaluations')) is int and result['evaluations'] in (
                PARENT_EVALUATIONS-1, PARENT_EVALUATIONS):
            return allow('declared_parent_limit', 'declared_parent_evaluation_limit')
        return deny('unknown_or_inconsistent', 'unadmitted_parent_status_or_evaluation_count')
    if arm == 'full':
        if status not in ('complete', 'cpu_limit'):
            return deny('algorithm_control_stop' if status in ('evaluation_limit', 'support_limit') else
                        'unknown_or_inconsistent', 'unadmitted_full_status')
        if type(stats.get('seed_evaluations')) is not int or stats['seed_evaluations'] != 4:
            return deny('insufficient_starts', 'full_control_requires_four_scored_starts')
        return allow('completed' if status == 'complete' else 'budget_stop', status)

    if not _order(parent_order) or len(parent_order) != len(result['order']):
        return deny('unknown_or_inconsistent', 'missing_or_invalid_frozen_parent_order')
    changed = list(result['order']) != list(parent_order)
    if not isinstance(ws, str) or any(value is not None and not isinstance(value, str) for value in (co, cr)):
        return deny('unknown_or_inconsistent', 'missing_or_malformed_window_status')
    if status not in ('complete', 'cpu_limit', 'evaluation_limit'):
        return deny('unknown_or_inconsistent', 'unadmitted_native_window_status')
    if co is not None and co not in CHILD_STOPS and co != 'complete':
        category = ('numerical_refusal' if co == 'child_numerical_limit' else
                    'ambiguous_support_refusal' if co == 'child_admission_or_support_limit' else
                    'runtime_unavailable' if co == 'owned_child_cpu_unavailable' else 'backend_failure')
        return deny(category, co)
    if cr is not None and (cr not in CHILD_REFUSALS or co != 'child_'+cr):
        return deny('numerical_refusal' if cr == 'numerical_limit' else
                    'ambiguous_support_refusal' if cr == 'admission_or_support_limit' else
                    'backend_failure', 'refused_or_inconsistent_child:'+cr)
    if ws in RUNTIME_FAILURES:
        return deny('runtime_unavailable', ws)
    if ws in REPLY_FAILURES or ws == 'child_launch_or_transport_failure':
        return deny('backend_failure', ws)
    if ws in DOMAIN_REFUSALS:
        return deny('domain_refusal', ws)
    if ws in ('numerical_limit', 'child_numerical_limit'):
        return deny('numerical_refusal', ws)
    if ws == 'child_admission_or_support_limit':
        return deny('ambiguous_support_refusal', ws)
    if ws in ('child_failure', 'child_worker_failure'):
        return deny('backend_failure', ws)
    if co is not None and (stats.get('backend_runtime_assets_verified') is not True or
            not _integer(stats.get('child_returncode'), -2**31) or
            (co == 'complete' and stats['child_returncode'] != 0)):
        return deny('unknown_or_inconsistent', 'unverified_assets_or_inconsistent_child_returncode')
    if (stats.get('parent_method') != 'magnitude' or type(stats.get('parent_evaluations')) is not int or
            stats['parent_evaluations'] != 1):
        return deny('unknown_or_inconsistent', 'supplied_parent_not_scored_once')

    def unchanged(category, reason):
        return deny('unverified_checkpoint', 'stop_selected_unverified_changed_order') if changed else allow(category, reason)

    if ws == 'support_limit':
        reason = stats.get('refusal_reason')
        if isinstance(reason, str) and reason in CAPACITY_REASONS:
            return unchanged('capacity_stop', reason)
        asset_error = isinstance(reason, str) and (reason.startswith('backend ') or
            reason.startswith('malformed backend ') or reason.startswith('installed Julia ') or
            reason == 'window backend requires Python 3.11+ tomllib')
        return deny('runtime_asset_failure' if asset_error else 'ambiguous_support_refusal',
                    reason if isinstance(reason, str) else 'missing_or_malformed_support_reason')
    if ws in PARENT_BUDGET_STOPS:
        return unchanged('budget_stop', ws)
    if ws in CHILD_STOPS:
        return unchanged('budget_stop', ws) if co == ws else deny('unknown_or_inconsistent', 'child_stop_outcome_differs')
    if co is None or not _checkpoint(stats):
        return deny('unverified_checkpoint', 'missing_or_invalid_parsed_checkpoint_metadata')
    if co == 'complete' and stats['child_checkpoint_status'] != 'complete':
        return deny('unknown_or_inconsistent', 'completed_child_has_incomplete_checkpoint')
    verified = (type(stats.get('native_verifications')) is int and stats['native_verifications'] == 1 and
                type(stats.get('native_verification_attempts')) is int and stats['native_verification_attempts'] == 1)
    if ws == 'complete' and co == 'complete' and not changed:
        return allow('completed', 'checked_complete_without_required_improvement')
    if ws == 'native_no_resolved_improvement' and not changed and verified:
        return allow('checked_no_improvement', ws)
    if changed and verified and stats['child_checked_updates'] >= 1:
        if ws == 'verified_improvement' and co == 'complete':
            return allow('verified_checkpoint', ws)
        if co in CHILD_STOPS and ws == 'verified_checkpoint_after_'+co:
            return allow('verified_checkpoint', ws)
    return deny('unknown_or_inconsistent', 'unknown_or_inconsistent_window_outcome')
