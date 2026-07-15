"""Dedicated single-qudit reset instruction.

``Reset`` re-initialises one qudit to ``|0⟩``: the qudit is measured in
the computational basis, the outcome is discarded, and the qudit is
prepared in ``|0⟩``.  As a quantum channel this is exactly the CPTP map
with Kraus operators ``{K_k = |0⟩⟨k|}_{k=0..d-1}``
(``Σ_k K_k† K_k = I``), so the class is implemented as a
:class:`KrausChannel` subclass:

* density-matrix backends (:class:`mqt.qudits.simulation.backends.DMSim`)
  execute it deterministically as ``ρ → Σ_k K_k ρ K_k†`` — i.e. trace
  out the qudit and re-initialise it to ``|0⟩⟨0|`` (no new backend code
  is required);
* state-vector backends that support mid-circuit stochastic channels
  (currently :class:`mqt.qudits.simulation.backends.TNSim`) execute it
  **stochastically** per run: the Born probabilities
  ``p_k = ⟨ψ|K_k†K_k|ψ⟩`` of the computational-basis outcomes of the
  target qudit are computed, an outcome ``k`` is sampled, the state is
  projected and renormalised, and the qudit is set to ``|0⟩``.  This is
  the quantum-trajectory (measure-and-discard) semantics, equivalent to
  tensorcircuit-ng's ``Circuit.general_kraus``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .kraus_channel import KrausChannel

if TYPE_CHECKING:
    from ..circuit import QuantumCircuit


class Reset(KrausChannel):
    """Reset one qudit to ``|0⟩`` (measure, discard outcome, re-prepare).

    Parameters
    ----------
    circuit:
        Parent quantum circuit.
    name:
        Instruction name (used for repr).
    target_qudit:
        Index of the qudit to reset.
    dimension:
        Dimension ``d`` of the target qudit.
    """

    def __init__(
        self,
        circuit: QuantumCircuit,
        name: str,
        target_qudit: int,
        dimension: int,
    ) -> None:
        kraus_ops = []
        for k in range(dimension):
            k_op = np.zeros((dimension, dimension), dtype=np.complex128)
            k_op[0, k] = 1.0
            kraus_ops.append(k_op)
        super().__init__(
            circuit=circuit,
            name=name,
            target_qudits=target_qudit,
            kraus_operators=kraus_ops,
            dimensions=dimension,
        )
        self.qasm_tag = "reset"
