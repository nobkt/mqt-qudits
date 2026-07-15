"""Kraus channel instruction for density-matrix backends.

This instruction represents a CPTP map ``ρ → Σ_k K_k ρ K_k†`` acting on a
specified subset of qudits.  It is **not** a unitary gate, so it has no
matrix representation; calling :meth:`to_matrix` raises an explicit
:class:`NotImplementedError`.  Only density-matrix backends that
recognise :class:`KrausChannel` (currently :class:`mqt.qudits.simulation.backends.DMSim`)
can execute circuits containing this instruction.

The CPTP (completely-positive trace-preserving) condition

    Σ_k K_k† K_k = I

is verified at construction time within an absolute tolerance.  No
heuristic correction is applied; an invalid set of Kraus operators
raises :class:`ValueError`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ...exceptions import CircuitError
from ..components.extensions.gate_types import GateTypes
from ..gate import Gate

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

    from ..circuit import QuantumCircuit


CPTP_TOLERANCE = 1e-9


class KrausChannel(Gate):
    """Non-unitary Kraus-channel instruction.

    Parameters
    ----------
    circuit:
        Parent quantum circuit.
    name:
        Instruction name (used for QASM/repr).
    target_qudits:
        Single qudit index or list of qudit indices the channel acts on.
    kraus_operators:
        List of Kraus matrices.  Each matrix has shape
        ``(D_local, D_local)`` with ``D_local = prod(dimensions)``.
        The convention for the leg ordering matches that of
        :class:`mqt.qudits.quantum_circuit.gates.CustomMulti`: row index
        ranges over the tensor product of the target qudits taken in
        the order given by ``target_qudits``.
    dimensions:
        Dimension of each target qudit (single ``int`` for a single
        qudit, list otherwise).
    """

    def __init__(
        self,
        circuit: QuantumCircuit,
        name: str,
        target_qudits: list[int] | int,
        kraus_operators: Sequence[NDArray[np.complex128]],
        dimensions: list[int] | int,
    ) -> None:
        # Determine gate type from number of target qudits (number of
        # target qudits, not d²-products).
        if isinstance(target_qudits, int):
            n_targets = 1
            local_dims: list[int] = [dimensions if isinstance(dimensions, int) else dimensions[0]]
        else:
            n_targets = len(target_qudits)
            if isinstance(dimensions, int):
                msg = "dimensions must be a list when target_qudits is a list"
                raise CircuitError(msg)
            local_dims = list(dimensions)

        if n_targets == 1:
            gate_type = GateTypes.SINGLE
        elif n_targets == 2:
            gate_type = GateTypes.TWO
        else:
            gate_type = GateTypes.MULTI

        d_local = int(np.prod(local_dims))

        # Validate and store Kraus operators (defensive copy).
        if len(kraus_operators) == 0:
            msg = "KrausChannel requires at least one Kraus operator"
            raise ValueError(msg)
        kraus_list: list[NDArray[np.complex128]] = []
        for idx, K in enumerate(kraus_operators):
            arr = np.asarray(K, dtype=np.complex128)
            if arr.shape != (d_local, d_local):
                msg = (
                    f"Kraus operator #{idx} has shape {arr.shape}; "
                    f"expected ({d_local},{d_local}) for dimensions {local_dims}."
                )
                raise ValueError(msg)
            kraus_list.append(arr.copy())

        # Verify CPTP: Σ K_k† K_k = I within tolerance.
        completeness = np.zeros((d_local, d_local), dtype=np.complex128)
        for K in kraus_list:
            completeness += K.conj().T @ K
        identity = np.eye(d_local, dtype=np.complex128)
        defect = np.linalg.norm(completeness - identity, ord="fro")
        if defect > CPTP_TOLERANCE:
            msg = (
                "Kraus operators are not trace-preserving: "
                f"||Σ K† K - I||_F = {defect:.3e} (tol {CPTP_TOLERANCE:.0e}). "
                "No heuristic correction is applied."
            )
            raise ValueError(msg)

        # Store kraus operators on the instance for backend retrieval.
        self._kraus_operators: list[NDArray[np.complex128]] = kraus_list
        self.is_kraus_channel: bool = True

        # Pass the first Kraus operator as `params` only for QASM bookkeeping;
        # this is **not** a unitary representation of the channel.
        super().__init__(
            circuit=circuit,
            name=name,
            gate_type=gate_type,
            target_qudits=target_qudits,
            dimensions=dimensions,
            control_set=None,
            params=kraus_list[0],
            qasm_tag="kraus",
        )

    @property
    def kraus_operators(self) -> list[NDArray[np.complex128]]:
        """Return the stored Kraus operators (no defensive copy)."""
        return self._kraus_operators

    @property
    def kraus_rank(self) -> int:
        return len(self._kraus_operators)

    def __array__(self) -> NDArray[np.complex128]:  # noqa: PLW3201
        msg = (
            "KrausChannel has no unitary matrix representation. "
            "Use a density-matrix backend (e.g. DMSim) instead."
        )
        raise NotImplementedError(msg)

    def to_matrix(self, identities: int = 0) -> NDArray[np.complex128]:  # noqa: ARG002
        msg = (
            "KrausChannel has no unitary matrix representation. "
            "Use a density-matrix backend (e.g. DMSim) instead."
        )
        raise NotImplementedError(msg)

    @staticmethod
    def validate_parameter(parameter: NDArray | None = None) -> bool:  # noqa: ARG004
        # Validation is performed in __init__ on the full kraus list.
        return True
