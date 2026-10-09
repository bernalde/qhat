"""
Pauli String Time Evolution Bloq

This module implements a Qualtran Bloq for time evolution under a single Pauli string Hamiltonian:
    U = exp(-i * coefficient * PauliString * t / ℏ)

where:
    - PauliString is a tensor product of Pauli operators (I, X, Y, Z)
    - coefficient is a scalar weight
    - t is the evolution time
    - ℏ is the reduced Planck constant (defaults to 1 in natural units)
"""

from functools import cached_property
from typing import TYPE_CHECKING, Dict

import attrs
import numpy as np
from scipy.linalg import expm

from cirq import DensePauliString as CirqDensePauliString
from cirq import GlobalPhaseGate, PauliStringPhasorGate
from qualtran import Bloq, BloqBuilder, QBit, Register, Side, Signature, SoquetT
from qualtran.cirq_interop import CirqGateAsBloq
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity

from qhat.common.pauli_utils import validate_pauli_string

if TYPE_CHECKING:
    import quimb.tensor as qtn


@attrs.frozen
class PauliStringEvolution(Bloq):
    """Time evolution operator for a single Pauli string Hamiltonian.

    Implements the unitary operator:
        U = exp(-i * coefficient * PauliString * time / hbar)

    Attributes:
        pauli_string: String of Pauli operators, e.g., "IXXYIIZ". Each character must be
            one of 'I', 'X', 'Y', or 'Z'.
        coefficient: Scalar coefficient multiplying the Pauli string in the Hamiltonian.
        time: Evolution time parameter.
        hbar: Reduced Planck constant (default 1.0 for natural units).

    Example:
        >>> pse = PauliStringEvolution(pauli_string="XYZ", coefficient=0.5, time=1.0)
        >>> U = pse.tensor_contract()  # Get the unitary matrix
        >>> from qualtran.cirq_interop.t_complexity_protocol import t_complexity
        >>> tc = t_complexity(pse)  # Get resource estimates
    """

    pauli_string: str
    coefficient: float
    time: float
    hbar: float = 1.0

    def __attrs_post_init__(self):
        """Validate the Pauli string contains only valid characters."""
        validate_pauli_string(self.pauli_string)

    @cached_property
    def signature(self) -> Signature:
        """Return the signature with a register of qubits matching the Pauli string length."""
        # Use shape=(n,) to create an array of n individual qubits, not a single n-bit register
        return Signature([Register('q', QBit(), shape=(self.num_qubits,), side=Side.THRU)])

    @cached_property
    def num_qubits(self) -> int:
        """Number of qubits the operator acts on."""
        return len(self.pauli_string)

    def _get_pauli_matrix(self) -> np.ndarray:
        """Get the matrix representation of the Pauli string."""
        # Build up the full Pauli matrix via tensor products
        pauli_matrices = {
            'I': np.array([[1, 0], [0, 1]], dtype=complex),
            'X': np.array([[0, 1], [1, 0]], dtype=complex),
            'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
            'Z': np.array([[1, 0], [0, -1]], dtype=complex),
        }

        result = pauli_matrices[self.pauli_string[0]]
        for char in self.pauli_string[1:]:
            result = np.kron(result, pauli_matrices[char])

        return result

    def _compute_unitary_matrix(self) -> np.ndarray:
        """Compute the unitary matrix for this Pauli string evolution."""
        # Build the Pauli string matrix
        P = self._get_pauli_matrix()

        # Compute the evolution operator: exp(-i * coefficient * P * time / hbar)
        angle = -1j * self.coefficient * self.time / self.hbar
        U = expm(angle * P)

        return U

    def my_tensors(self, incoming, outgoing):
        """Supply the exact unitary through Qualtran's tensor protocol.

        Each shaped QBit connection has one bit index, zero. Preserve the global
        phase so tensor simulation and controlled circuit evolution agree.
        """
        import quimb.tensor as qtn

        return [qtn.Tensor(
            data=self._compute_unitary_matrix().reshape([2] * (2 * self.num_qubits)),
            inds=[(outgoing['q'][i], 0) for i in range(self.num_qubits)]
                 + [(incoming['q'][i], 0) for i in range(self.num_qubits)],
            tags=[str(self)],
        )]

    def _build_cirq_gate(self):
        """Build the corresponding Cirq gate for resource estimation."""
        # The parity gate acts only on the nonidentity support. Passing idle
        # positions through a parity decomposition can rotate the wrong Pauli.
        support_string = self.pauli_string.replace('I', '')
        if not support_string:
            return GlobalPhaseGate(np.exp(-1j * self.coefficient * self.time / self.hbar))
        cirq_pauli = CirqDensePauliString(support_string, coefficient=1.0)

        # Cirq phases the -1 and +1 eigenspaces by exp(i*pi*exponent).
        # Setting both exponents implements exp(-i*theta*P) exactly. Using the
        # default positive-eigenspace exponent of zero would introduce an
        # extra exp(i*theta), which changes controlled evolution and energies.
        half_turns = self.coefficient * self.time / (self.hbar * np.pi)
        return PauliStringPhasorGate(
            cirq_pauli, exponent_neg=half_turns, exponent_pos=-half_turns,
        )

    def build_composite_bloq(self, bb: BloqBuilder, **soqs: SoquetT) -> Dict[str, SoquetT]:
        """Decompose into a Cirq gate for resource estimation.

        This decomposition is needed for pyLIQTR and controlled versions to work properly.
        """
        # Get the Cirq gate representation
        cirq_gate = self._build_cirq_gate()

        # Wrap it as a Bloq
        gate_bloq = CirqGateAsBloq(cirq_gate)

        # Keep identity positions as idle wires. CirqGateAsBloq uses a scalar
        # QBit for support size one and a shaped register for larger support.
        support = [i for i, char in enumerate(self.pauli_string) if char != 'I']
        q = np.array(soqs['q'], dtype=object, copy=True)
        if not support:
            bb.add(gate_bloq)
        elif len(support) == 1:
            q[support[0]] = bb.add(gate_bloq, q=q[support[0]])
        else:
            q[support] = bb.add(gate_bloq, q=q[support])
        return {'q': q}

    def _t_complexity_(self) -> TComplexity:
        """Return the T-complexity for resource estimation."""
        # Build the internal Cirq gate representation and get its complexity
        cirq_gate = self._build_cirq_gate()
        bloq = CirqGateAsBloq(cirq_gate)
        return t_complexity(bloq)

    def __str__(self) -> str:
        return f"exp(-i * {self.coefficient} * {self.pauli_string} * {self.time} / {self.hbar})"

    def __repr__(self) -> str:
        return (f"PauliStringEvolution(pauli_string='{self.pauli_string}', "
                f"coefficient={self.coefficient}, time={self.time}, hbar={self.hbar})")
