from __future__ import annotations

import operator
from functools import reduce
from typing import TYPE_CHECKING

import numpy as np
import tensornetwork as tn  # type: ignore[import-not-found]
from typing_extensions import Unpack

from ...quantum_circuit.components.extensions.gate_types import GateTypes
from ...quantum_circuit.gates.kraus_channel import KrausChannel
from ..jobs import Job, JobResult
from .backendv2 import Backend
from .stochastic_sim import stochastic_simulation

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

    from ...quantum_circuit import QuantumCircuit
    from ...quantum_circuit.gate import Gate
    from .. import MQTQuditProvider
    from ..noise_tools import NoiseModel


class TNSim(Backend):
    def __init__(
        self,
        provider: MQTQuditProvider,
        name: str | None = None,
        description: str | None = None,
        **fields: Unpack[Backend.DefaultOptions],
    ) -> None:
        super().__init__(provider, name=name, description=description, **fields)

    def __noise_model(self) -> NoiseModel | None:
        return self.noise_model

    def run(self, circuit: QuantumCircuit, **options: Unpack[Backend.DefaultOptions]) -> Job:
        job = Job(self)
        self._options.update(options)
        self.noise_model: NoiseModel | None = self._options.get("noise_model", None)
        self.shots = self._options.get("shots", 50)
        self.memory = self._options.get("memory", False)
        self.full_state_memory = self._options.get("full_state_memory", False)
        self.file_path = self._options.get("file_path", None)
        self.file_name = self._options.get("file_name", None)
        self._rng = np.random.default_rng(self._options.get("seed", None))

        if self.noise_model is not None:
            assert self.shots >= 50, "Number of shots should be above 50"
            job.set_result(JobResult(state_vector=self.execute(circuit), counts=stochastic_simulation(self, circuit)))
        else:
            job.set_result(JobResult(state_vector=self.execute(circuit), counts=[]))

        return job

    def execute(self, circuit: QuantumCircuit, noise_model: NoiseModel | None = None) -> NDArray[np.complex128]:  # noqa: ARG002
        self.system_sizes = circuit.dimensions
        self.circ_operations = circuit.instructions

        if any(isinstance(op, KrausChannel) for op in self.circ_operations):
            # Mid-circuit non-unitary channels (KrausChannel / Reset):
            # stochastic single-trajectory execution (measure-and-discard
            # semantics).  Each `execute` call yields ONE trajectory; the
            # returned state vector is a sample, not an average.
            psi_t = self.__evolve_with_channels(self.system_sizes, self.circ_operations)
            state_size = reduce(operator.mul, self.system_sizes, 1)
            return psi_t.reshape(1, state_size)

        result = self.__contract_circuit(self.system_sizes, self.circ_operations)

        result = np.transpose(result.tensor, list(range(len(self.system_sizes))))

        state_size = reduce(operator.mul, self.system_sizes, 1)
        return result.reshape(1, state_size)

    @staticmethod
    def _apply_local_matrix_to_state(
        psi_t: NDArray[np.complex128],
        op_matrix: NDArray[np.complex128],
        qudits: Sequence[int],
        dims: Sequence[int],
    ) -> NDArray[np.complex128]:
        """Apply a ``(D_local, D_local)`` matrix to the given qudit legs of ``psi_t``.

        ``psi_t`` is the state reshaped to ``(*dims,)``.  The matrix acts on
        the tensor product of the target qudits in the order given by
        ``qudits`` (same convention as :class:`KrausChannel` /
        :class:`CustomMulti`).
        """
        n = len(dims)
        k = len(qudits)
        local_dims = [dims[q] for q in qudits]
        op_t = op_matrix.reshape(*local_dims, *local_dims)
        out = np.tensordot(op_t, psi_t, axes=(list(range(k, 2 * k)), list(qudits)))
        remaining = [a for a in range(n) if a not in qudits]
        label_to_pos: dict[int, int] = {}
        for new_pos, label in enumerate(list(qudits) + remaining):
            label_to_pos[label] = new_pos
        perm = [label_to_pos[i] for i in range(n)]
        return np.transpose(out, perm)

    def _apply_kraus_stochastic(
        self,
        psi_t: NDArray[np.complex128],
        channel: KrausChannel,
        dims: Sequence[int],
    ) -> NDArray[np.complex128]:
        """Sample one Kraus branch with Born probability and project.

        ``p_k = ⟨ψ|K_k†K_k|ψ⟩``; the state collapses to
        ``K_k|ψ⟩/√p_k``.  This realises the measure-and-discard
        (quantum-trajectory) semantics of a CPTP channel on a pure state,
        equivalent to tensorcircuit-ng's ``Circuit.general_kraus``.
        """
        target = channel.target_qudits
        qudits: tuple[int, ...] = (target,) if isinstance(target, int) else tuple(target)
        branches: list[NDArray[np.complex128]] = []
        probs: list[float] = []
        for k_op in channel.kraus_operators:
            cand = self._apply_local_matrix_to_state(psi_t, k_op, qudits, dims)
            branches.append(cand)
            probs.append(float(np.real(np.vdot(cand, cand))))
        prob_arr = np.array(probs, dtype=np.float64)
        total = float(prob_arr.sum())
        # The channel is CPTP (verified at construction), so on a normalised
        # state the probabilities must sum to 1 up to round-off.
        if not np.isclose(total, 1.0, atol=1e-8):
            msg = (
                "Kraus branch probabilities do not sum to 1 "
                f"(got {total:.6e}); the input state may not be normalised."
            )
            raise ValueError(msg)
        rng = getattr(self, "_rng", None)
        if rng is None:
            rng = np.random.default_rng()
            self._rng = rng
        idx = int(rng.choice(len(prob_arr), p=prob_arr / total))
        p_sel = probs[idx]
        return branches[idx] / np.sqrt(p_sel)

    def __evolve_with_channels(
        self, system_sizes: list[int], operations: Sequence[Gate]
    ) -> NDArray[np.complex128]:
        """Single stochastic trajectory through a circuit containing channels.

        Unitary sub-sequences between channels are contracted with the
        regular tensor-network path (starting from the current state);
        each :class:`KrausChannel` (including :class:`Reset`) is applied
        stochastically via Born sampling.
        """
        psi_t: NDArray[np.complex128] | None = None  # None = |0…0⟩ product state
        pending: list[Gate] = []

        def flush(state: NDArray[np.complex128] | None) -> NDArray[np.complex128] | None:
            if not pending:
                return state
            result = self.__contract_circuit(system_sizes, list(pending), initial_state=state)
            pending.clear()
            return np.transpose(result.tensor, list(range(len(system_sizes))))

        for op in operations:
            if isinstance(op, KrausChannel):
                psi_t = flush(psi_t)
                if psi_t is None:
                    psi_t = np.zeros(tuple(system_sizes), dtype=np.complex128)
                    psi_t[(0,) * len(system_sizes)] = 1.0
                psi_t = self._apply_kraus_stochastic(psi_t, op, system_sizes)
            else:
                pending.append(op)
        psi_t = flush(psi_t)
        if psi_t is None:
            psi_t = np.zeros(tuple(system_sizes), dtype=np.complex128)
            psi_t[(0,) * len(system_sizes)] = 1.0
        return psi_t

    @staticmethod
    def __apply_gate(qudit_edges: tn.Edge, gate: NDArray, operating_qudits: list[int]) -> None:
        op = tn.Node(gate)
        for i, bit in enumerate(operating_qudits):
            tn.connect(qudit_edges[bit], op[i])
            qudit_edges[bit] = op[i + len(operating_qudits)]

    def __contract_circuit(
        self,
        system_sizes: list[int],
        operations: Sequence[Gate],
        initial_state: NDArray[np.complex128] | None = None,
    ) -> tn.network_components.AbstractNode:
        all_nodes: Sequence[tn.network_components.AbstractNode] = []

        with tn.NodeCollection(all_nodes):
            if initial_state is None:
                state_nodes = []
                for s in system_sizes:
                    z = [0] * s
                    z[0] = 1
                    state_nodes.append(tn.Node(np.array(z, dtype="complex")))

                qudits_legs = [node[0] for node in state_nodes]
            else:
                init_node = tn.Node(np.asarray(initial_state, dtype="complex"))
                qudits_legs = [init_node[i] for i in range(len(system_sizes))]

            for op in operations:
                op_matrix = op.to_matrix(identities=1)
                lines = op.reference_lines

                if op.gate_type == GateTypes.SINGLE:
                    op_matrix = op_matrix.T
                    # op_matrix = op_matrix.reshape((system_sizes[lines[0]], system_sizes[lines[0]]))

                elif op.gate_type == GateTypes.TWO and not op.is_long_range:
                    op_matrix = op_matrix.T
                    lines = lines.copy()
                    lines.sort()

                    op_matrix = op_matrix.reshape((
                        system_sizes[lines[0]],
                        system_sizes[lines[1]],
                        system_sizes[lines[0]],
                        system_sizes[lines[1]],
                    ))

                elif op.is_long_range or op.gate_type == GateTypes.MULTI:
                    op_matrix = op_matrix.T
                    minimum_line, maximum_line = min(lines), max(lines)
                    interested_lines = list(range(minimum_line, maximum_line + 1))
                    inputs_outputs_legs = [system_sizes[i] for i in interested_lines] * 2
                    op_matrix = op_matrix.reshape(tuple(inputs_outputs_legs))
                    lines = interested_lines

                self.__apply_gate(qudits_legs, op_matrix, lines)

        return tn.contractors.auto(all_nodes, output_edge_order=qudits_legs)
