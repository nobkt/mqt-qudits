"""Density-matrix simulator backend.

This backend evolves a density matrix ``ρ`` directly under the
instructions of a :class:`mqt.qudits.quantum_circuit.QuantumCircuit`.
Unlike :class:`TNSim` and :class:`MISim` (which contract a pure state
vector), this backend stores ``ρ`` as a ``(D, D)`` complex matrix where
``D = prod(dimensions)``.  This is mathematically equivalent to
running an open-system quantum simulation directly: it is the standard
"density matrix simulator" concept (e.g. ``DensityMatrixSimulator`` in
Qiskit Aer / Cirq).

Two kinds of instructions are supported:

* **Unitary gates** (``CustomOne``, ``CustomTwo``, ``CustomMulti`` and
  any other :class:`Gate` subclass that returns a unitary matrix from
  ``to_matrix(identities=1)``): applied as ``ρ → U_S ρ U_S†`` where
  ``U_S`` is the unitary embedded on the target qudit subset ``S``.

* **Kraus channels** (:class:`KrausChannel`): applied as
  ``ρ → Σ_k K_k ρ K_k†``.  A density-matrix backend is required because
  a single CPTP channel cannot in general be realised by a fixed
  unitary acting on the system alone.

Initial state
-------------
By default ``ρ_0 = |0…0⟩⟨0…0|``.  An arbitrary positive
trace-1 Hermitian initial density matrix can be supplied via the
backend option ``initial_density_matrix``.

Capacity
--------
For ``D = prod(dimensions)`` the storage is ``16 D²`` bytes.  For the
Scenario 5 GKSL system (``N=4`` qutrits, ``D = 81``) this is
``≈ 100 KB``: trivially tractable.  Per-instruction cost is
``O(D · D_local²)`` for unitary gates and ``O(rank · D · D_local²)`` for
Kraus channels, where ``D_local = prod(d_q for q in target_qudits)``.

This backend is *not* a substitute for the state-vector backends
``tnsim`` / ``misim``: it simulates the density matrix directly and
does not produce a state-vector readout.  Its ``JobResult.state_vector``
attribute holds the **flattened density matrix** for compatibility
with the existing :class:`JobResult` API; the unflattened ``ρ`` is also
exposed via :attr:`JobResult.density_matrix`.
"""

from __future__ import annotations

import operator
from functools import reduce
from typing import TYPE_CHECKING

import numpy as np
from typing_extensions import Unpack

from ...quantum_circuit.gates.kraus_channel import KrausChannel
from ..jobs import Job, JobResult
from .backendv2 import Backend

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

    from ...quantum_circuit import QuantumCircuit
    from ...quantum_circuit.gate import Gate
    from .. import MQTQuditProvider
    from ..noise_tools import NoiseModel


# -----------------------------------------------------------------------
# Tensor helpers
# -----------------------------------------------------------------------


def _apply_local_op_one_side(
    rho_tensor: NDArray[np.complex128],
    op_matrix: NDArray[np.complex128],
    qudits: Sequence[int],
    dims: Sequence[int],
    side: str,
) -> NDArray[np.complex128]:
    """Multiply rho on row (side='row') or column (side='col') by op_matrix.

    Parameters
    ----------
    rho_tensor:
        Density matrix reshaped to ``(*dims, *dims)`` — first ``N`` axes
        are row legs, last ``N`` axes are column legs.
    op_matrix:
        ``(D_local, D_local)`` matrix where
        ``D_local = prod(dims[q] for q in qudits)``.  For the row pass we
        compute ``ρ_new = (op ⊗ I) ρ`` (acts on row legs); for the column
        pass we compute ``ρ_new = ρ (op ⊗ I)`` (acts on column legs) — to
        realise ``ρ → A ρ A†`` one therefore calls this with ``op=A`` on
        side ``'row'`` and then with ``op=A.conj()`` on side ``'col'``.
    qudits:
        Tuple of qudit indices the operator acts on, in the order
        matching the leg ordering of ``op_matrix``.
    dims:
        Full list of qudit dimensions.
    side:
        ``'row'`` or ``'col'``.
    """
    n = len(dims)
    k = len(qudits)
    local_dims = [dims[q] for q in qudits]

    op_t = op_matrix.reshape(*local_dims, *local_dims)

    if side == "row":
        contract_axes = list(qudits)
    elif side == "col":
        contract_axes = [n + q for q in qudits]
    else:
        msg = f"side must be 'row' or 'col', got {side!r}"
        raise ValueError(msg)

    # tensordot contracts op_t's "in" legs (positions k..2k-1) with rho_t
    # legs at positions `contract_axes`.  The result has axes in the
    # order: (op_out_0, ..., op_out_{k-1}, *remaining_rho_axes_in_order).
    out = np.tensordot(op_t, rho_tensor, axes=(list(range(k, 2 * k)), contract_axes))

    remaining = [a for a in range(2 * n) if a not in contract_axes]
    # Build inverse permutation: we want axis labels [0, 1, ..., 2N-1] in order.
    label_to_pos: dict[int, int] = {}
    for new_pos, label in enumerate(list(contract_axes) + remaining):
        label_to_pos[label] = new_pos
    perm = [label_to_pos[i] for i in range(2 * n)]
    return np.transpose(out, perm)


def apply_unitary_to_density(
    rho: NDArray[np.complex128],
    u_matrix: NDArray[np.complex128],
    qudits: Sequence[int],
    dims: Sequence[int],
) -> NDArray[np.complex128]:
    """Return ``U_S ρ U_S†`` where ``U_S`` is ``u_matrix`` on the target qudits."""
    n = len(dims)
    d_total = int(reduce(operator.mul, dims, 1))
    rho_t = rho.reshape(*([d for d in dims] * 2)) if n > 0 else rho
    rho_t = _apply_local_op_one_side(rho_t, u_matrix, qudits, dims, "row")
    rho_t = _apply_local_op_one_side(rho_t, u_matrix.conj(), qudits, dims, "col")
    return rho_t.reshape(d_total, d_total)


def apply_kraus_to_density(
    rho: NDArray[np.complex128],
    kraus_ops: Sequence[NDArray[np.complex128]],
    qudits: Sequence[int],
    dims: Sequence[int],
) -> NDArray[np.complex128]:
    """Return ``Σ_k K_k ρ K_k†`` for Kraus ops acting on the target qudits."""
    n = len(dims)
    d_total = int(reduce(operator.mul, dims, 1))
    accumulator = np.zeros((d_total, d_total), dtype=np.complex128)
    for k_op in kraus_ops:
        rho_t = rho.reshape(*([d for d in dims] * 2)) if n > 0 else rho
        rho_t = _apply_local_op_one_side(rho_t, k_op, qudits, dims, "row")
        rho_t = _apply_local_op_one_side(rho_t, k_op.conj(), qudits, dims, "col")
        accumulator += rho_t.reshape(d_total, d_total)
    return accumulator


# -----------------------------------------------------------------------
# JobResult extension
# -----------------------------------------------------------------------


class DensityMatrixJobResult(JobResult):
    """JobResult that additionally exposes the final density matrix."""

    def __init__(self, density_matrix: NDArray[np.complex128]) -> None:
        flat = density_matrix.reshape(1, -1)
        super().__init__(state_vector=flat, counts=[])
        self._density_matrix = density_matrix

    @property
    def density_matrix(self) -> NDArray[np.complex128]:
        return self._density_matrix

    def get_density_matrix(self) -> NDArray[np.complex128]:
        return self._density_matrix


# -----------------------------------------------------------------------
# DMSim backend
# -----------------------------------------------------------------------


class DMSim(Backend):
    """Density-matrix simulator backend.

    Recognised options (in addition to :class:`Backend.DefaultOptions`):

    * ``initial_density_matrix``: ``(D, D)`` complex array used as
      ``ρ_0`` (default: ``|0…0⟩⟨0…0|``).  No CPTP / Hermiticity check
      is performed beyond shape; the caller is responsible for supplying
      a physical density matrix.

    The backend ignores any ``noise_model`` option (raises if one is
    supplied — noise must be supplied as :class:`KrausChannel`
    instructions so that the simulator and the user agree on what the
    physics is).
    """

    def __init__(
        self,
        provider: MQTQuditProvider,
        name: str | None = None,
        description: str | None = None,
        **fields: Unpack[Backend.DefaultOptions],
    ) -> None:
        super().__init__(provider, name=name, description=description, **fields)

    def run(self, circuit: QuantumCircuit, **options: Unpack[Backend.DefaultOptions]) -> Job:
        job = Job(self)
        self._options.update(options)
        if self._options.get("noise_model", None) is not None:
            msg = (
                "DMSim does not accept a NoiseModel option. "
                "Add noise as KrausChannel instructions in the circuit instead."
            )
            raise ValueError(msg)
        rho = self.execute(circuit, initial_density_matrix=options.get("initial_density_matrix"))
        job.set_result(DensityMatrixJobResult(rho))
        return job

    def execute(
        self,
        circuit: QuantumCircuit,
        noise_model: NoiseModel | None = None,  # noqa: ARG002
        initial_density_matrix: NDArray[np.complex128] | None = None,
    ) -> NDArray[np.complex128]:
        dims: list[int] = list(circuit.dimensions)
        n = len(dims)
        d_total = int(reduce(operator.mul, dims, 1))

        if initial_density_matrix is None:
            psi0 = np.zeros(d_total, dtype=np.complex128)
            psi0[0] = 1.0
            rho = np.outer(psi0, psi0.conj())
        else:
            rho = np.asarray(initial_density_matrix, dtype=np.complex128)
            if rho.shape != (d_total, d_total):
                msg = (
                    f"initial_density_matrix has shape {rho.shape}; "
                    f"expected ({d_total},{d_total}) for circuit dimensions "
                    f"{dims}."
                )
                raise ValueError(msg)
            rho = rho.copy()

        for instruction in circuit.instructions:
            rho = self._apply_instruction(rho, instruction, dims, n)

        return rho

    @staticmethod
    def _apply_instruction(
        rho: NDArray[np.complex128],
        instruction: Gate,
        dims: list[int],
        n_qudits: int,  # noqa: ARG004
    ) -> NDArray[np.complex128]:
        target = instruction.target_qudits
        if isinstance(target, int):
            qudits: tuple[int, ...] = (target,)
        else:
            qudits = tuple(target)

        if isinstance(instruction, KrausChannel):
            return apply_kraus_to_density(rho, instruction.kraus_operators, qudits, dims)

        # Unitary gate path.
        u_matrix = instruction.to_matrix(identities=0)

        # Some Gate subclasses return the matrix in transposed form
        # depending on the leg convention.  Empirically (matching the
        # convention used by TNSim — see tnsim.py: op_matrix = op_matrix.T
        # for SINGLE / TWO before applying to legs), the canonical
        # action of `op_matrix` on the row legs in the order given by
        # `target_qudits` is achieved without an additional transpose
        # because TNSim transposes back into "legs out, legs in" form.
        # For DMSim we use the matrix directly as ``out, in`` operator
        # acting on the tensor product of target qudits in the order
        # ``target_qudits`` (the same convention KrausChannel uses).
        return apply_unitary_to_density(rho, u_matrix, qudits, dims)
