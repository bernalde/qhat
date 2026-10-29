"""Plan a chronological Strang formula under a Pauli-rotation budget.

Each unit is a nonempty group of mutually commuting Pauli terms. The caller
must establish that condition: this module sees group sizes, not operators.
Units follow the supplied chronological order; the first is outermost and the
last is central. Adjacent copies of the first unit are fused across steps.
Only these declared boundaries are fused; no further commutation optimization
is attempted. Budget counts emitted nonidentity Pauli rotations, not depth.

This maintained planner does not modify archived benchmark drivers or prove
the Python implementation correct in Lean.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from operator import index
from typing import Iterable, Iterator


def _integer(value: object, name: str, *, minimum: int) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an integer (not bool)")
    try:
        result = index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer (not bool)") from error
    if result < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return result


def _lengths(unit_lengths: Iterable[int]) -> tuple[int, ...]:
    lengths = tuple(
        _integer(length, f"unit_lengths[{index}]", minimum=1)
        for index, length in enumerate(unit_lengths)
    )
    if not lengths:
        raise ValueError("unit_lengths must contain at least one nonempty unit")
    return lengths


@dataclass(frozen=True)
class StrangPlan:
    """Immutable budget decision; use :meth:`iter_pulses` for chronological factors.

    A pulse ``(term_index, numerator)`` means the rotation for that flattened
    term has duration ``T * numerator / duration_denominator``. An infeasible
    plan has zero steps, zero rotations, denominator None, and no pulses; it
    must not be propagated. A singleton plan uses one exact commuting-group
    exponential, even if the budget could pay for redundant repetitions.
    """

    unit_lengths: tuple[int, ...]
    budget: int
    steps: int
    rotations: int

    @property
    def feasible(self) -> bool:
        return self.steps > 0

    @property
    def duration_denominator(self) -> int | None:
        return 2 * self.steps if self.feasible else None

    def iter_pulses(self) -> Iterator[tuple[int, int]]:
        """Yield fused pulses, or nothing when the budget is infeasible."""
        if not self.feasible:
            return iter(())
        return iter_strang_pulses(self.unit_lengths, self.steps)


def iter_strang_pulses(
    unit_lengths: Iterable[int], steps: int
) -> Iterator[tuple[int, int]]:
    """Return integer-duration pulses for a positive number of Strang steps.

    Validation happens at the call, before iteration. Each term receives
    total numerator ``2 * steps``. For one commuting unit, emit each term
    once with that numerator. For multiple units, preserve within-unit term
    order on both passes, which is valid under the commutation precondition.
    ``steps=0`` is invalid; handle infeasibility through ``StrangPlan.iter_pulses``.
    """
    lengths = _lengths(unit_lengths)
    count = _integer(steps, "steps", minimum=1)
    offsets = [0]
    for length in lengths:
        offsets.append(offsets[-1] + length)

    def emit() -> Iterator[tuple[int, int]]:
        if len(lengths) == 1:
            for term in range(lengths[0]):
                yield term, 2 * count
            return

        def unit(index: int, numerator: int) -> Iterator[tuple[int, int]]:
            for term in range(offsets[index], offsets[index + 1]):
                yield term, numerator

        yield from unit(0, 1)
        for step in range(count):
            for index in range(1, len(lengths) - 1):
                yield from unit(index, 1)
            yield from unit(len(lengths) - 1, 2)
            for index in range(len(lengths) - 2, 0, -1):
                yield from unit(index, 1)
            yield from unit(0, 1 if step == count - 1 else 2)

    return emit()


def plan_strang_schedule(
    unit_lengths: Iterable[int], budget: int
) -> StrangPlan:
    """Choose the largest feasible step count, with the commuting exception.

    Inputs are nonempty positive integer sizes and a nonnegative integer
    budget. Booleans and fractional numbers are rejected. With M >= 2 and
    U = sum(unit_lengths), the fused cost is

        steps * (2*U - unit_lengths[0] - unit_lengths[-1]) + unit_lengths[0].

    A single unit instead costs U rotations independently of steps: select
    one step if affordable. No positive step fitting means an infeasible
    plan with no pulses. Empty targets are deliberately outside this API.
    """
    lengths = _lengths(unit_lengths)
    available = _integer(budget, "budget", minimum=0)
    terms = sum(lengths)
    if len(lengths) == 1:
        steps = 1 if available >= terms else 0
        rotations = terms if steps else 0
    else:
        slope = 2 * terms - lengths[0] - lengths[-1]
        steps = max(0, (available - lengths[0]) // slope)
        rotations = steps * slope + lengths[0] if steps else 0
    return StrangPlan(lengths, available, steps, rotations)


def main(argv: list[str] | None = None) -> int:
    """Print a standalone, reviewable plan; this CLI does not propagate states."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-lengths", type=int, nargs="+", required=True)
    parser.add_argument("--budget", type=int, required=True)
    args = parser.parse_args(argv)
    try:
        plan = plan_strang_schedule(args.unit_lengths, args.budget)
    except (TypeError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps({
        "unit_lengths": plan.unit_lengths,
        "rotation_budget": plan.budget,
        "feasible": plan.feasible,
        "steps": plan.steps,
        "pauli_rotations": plan.rotations,
        "duration_denominator": plan.duration_denominator,
        "pulses": list(plan.iter_pulses()),
        "orientation": "first unit outermost; last unit central",
        "pulse_semantics": "[flattened term index, duration numerator]; duration = T * numerator / denominator",
        "precondition": "each declared unit consists of mutually commuting Pauli terms",
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
