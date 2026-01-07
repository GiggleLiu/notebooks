#set text(font: "New Computer Modern", size: 11pt)
#set page(margin: 1in)
#set par(justify: true)
#set heading(numbering: "1.")
#import "@preview/quill:0.6.0": *

#align(center)[
  #text(size: 18pt, weight: "bold")[Variational Quantum Compilation Challenge]
  #v(0.5em)
  #line(length: 50%, stroke: 0.5pt)
]

#v(1em)

= Background

Variational quantum compilation (VQC) approximates a target unitary $U_"target"$ with a parameterized circuit $V(bold(theta))$ by maximizing the fidelity #cite(<khatri2019qaqc>):
$ F = |"Tr"(U_"target"^dagger V(bold(theta)))|^2 \/ d^2 $

The compilation problem involves both *parameter optimization* (continuous, tractable) and *structure optimization* (discrete, NP-hard) #cite(<botea2018complexity>).

= Shor's Algorithm: The Bottleneck

Shor's algorithm for factoring $N$ consists of:

#table(
  columns: (auto, auto, auto),
  inset: 6pt,
  align: (left, center, left),
  [*Component*], [*Complexity*], [*Notes*],
  [Hadamard layer], [$O(n)$], [Initial superposition],
  [*Modular exponentiation*], [*$O(n^3)$*], [*Bottleneck: $2n$ controlled mod-mult*],
  [Inverse QFT], [$O(n^2)$], [$n(n-1)/2$ controlled rotations],
)

The modular exponentiation computes $a^x mod N$ via $2n$ controlled modular multiplications, each requiring $O(n^2)$ gates #cite(<beauregard2003>). This is the dominant cost.

= Challenge

== Target Unitary

Compile the *4-qubit modular multiplication* from Shor's algorithm for factoring $N = 15$:
$ U_7: |y chevron.r arrow.r |7 dot y mod 15 chevron.r $
with $a = 7$, $N = 15$ (the smallest non-trivial Shor instance #cite(<vandersypen2001>)).

This is a *permutation matrix* with period $r = 4$ (since $7^4 equiv 1 mod 15$).

== Circuit Structure

The target unitary $U_7$ can be decomposed into layers of SU(4) gates:

#figure(
  quantum-circuit(
    lstick($|q_0 chevron.r$), 1, mqgate($G_1$, n: 2), 1, 1, mqgate($G_4$, n: 2), 1, rstick($|q'_0 chevron.r$), [\ ],
    lstick($|q_1 chevron.r$), 1, 1, mqgate($G_2$, n: 2), mqgate($G_3$, n: 2), 1, 1, rstick($|q'_1 chevron.r$), [\ ],
    lstick($|q_2 chevron.r$), 1, 1, 1, 1, 1, 1, rstick($|q'_2 chevron.r$), [\ ],
    lstick($|q_3 chevron.r$), 1, 1, 1, 1, 1, 1, rstick($|q'_3 chevron.r$),
  ),
  caption: [Example ansatz structure with 4 SU(4) gates. Each $G_i in "SU"(4)$ acts on adjacent qubits.]
)

*Connectivity*: Linear chain `q0 ── q1 ── q2 ── q3`

== Objective

#align(center)[
  #rect(inset: 10pt, stroke: 1.5pt + blue, radius: 3pt)[
    Find circuit $V(bold(theta))$ with *minimum SU(4) count* achieving $F >= 0.99$
  ]
]

= Target Specification

== Permutation Table

The modular multiplication $|y chevron.r arrow.r |7 y mod 15 chevron.r$ maps:

#table(
  columns: 8,
  inset: 4pt,
  align: center,
  [$y$], [0], [1], [2], [3], [4], [5], [6], [7],
  [$7y mod 15$], [0], [7], [14], [6], [13], [5], [12], [4],
  [$y$], [8], [9], [10], [11], [12], [13], [14], [15],
  [$7y mod 15$], [11], [3], [10], [2], [9], [1], [8], [15],
)

*Key cycle*: $1 arrow.r 7 arrow.r 4 arrow.r 13 arrow.r 1$ (period 4)

== Target Unitary (Qiskit)

#text(size: 9pt)[
```python
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, process_fidelity
from qiskit.circuit.library import UnitaryGate

def build_mod_mult_unitary(a: int, N: int, n_qubits: int) -> np.ndarray:
    """Build the modular multiplication unitary U_a: |y⟩ → |a*y mod N⟩"""
    dim = 2**n_qubits
    U = np.zeros((dim, dim), dtype=complex)
    for y in range(dim):
        new_y = (a * y) % N if y < N else y
        U[new_y, y] = 1
    return U

# Build target unitary
U_target = build_mod_mult_unitary(a=7, N=15, n_qubits=4)
target_op = Operator(U_target)

# Verify properties
print(f"Target is unitary: {target_op.is_unitary()}")
print(f"Period is 4: {np.allclose(np.linalg.matrix_power(U_target, 4), np.eye(16))}")
```
]

= Verification Program (Qiskit)

#text(size: 9pt)[
```python
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, process_fidelity
from qiskit.circuit.library import UnitaryGate
import numpy as np

def verify_circuit(qc: QuantumCircuit, U_target: np.ndarray) -> dict:
    """
    Verify a Qiskit circuit against the target unitary.
    
    Args:
        qc: Qiskit QuantumCircuit with the compiled circuit
        U_target: Target unitary as numpy array
    
    Returns:
        dict with fidelity, gate_count, is_valid
    """
    # Get circuit unitary
    circuit_op = Operator(qc)
    target_op = Operator(U_target)
    
    # Compute process fidelity
    fidelity = process_fidelity(circuit_op, target_op)
    
    # Count 2-qubit gates
    gate_counts = qc.count_ops()
    two_qubit_gates = sum(v for k, v in gate_counts.items() 
                         if k in ['cx', 'cz', 'swap', 'unitary'] 
                         and qc.data[0].operation.num_qubits == 2)
    
    # Count custom unitaries (SU4 gates)
    su4_count = sum(1 for inst in qc.data 
                    if isinstance(inst.operation, UnitaryGate) 
                    and inst.operation.num_qubits == 2)
    
    return {
        "fidelity": fidelity,
        "su4_count": su4_count,
        "total_2q_gates": two_qubit_gates,
        "is_valid": fidelity >= 0.99
    }

def build_ansatz(n_layers: int = 2) -> QuantumCircuit:
    """Build a brick-wall ansatz with SU(4) gates."""
    qc = QuantumCircuit(4)
    
    for layer in range(n_layers):
        # Even layer: gates on (0,1), (2,3)
        # Odd layer: gates on (1,2)
        if layer % 2 == 0:
            for q in [0, 2]:
                if q + 1 < 4:
                    # Placeholder for SU(4) - use random unitary
                    from scipy.stats import unitary_group
                    G = unitary_group.rvs(4)
                    qc.unitary(G, [q, q+1], label=f'G{layer}_{q}')
        else:
            from scipy.stats import unitary_group
            G = unitary_group.rvs(4)
            qc.unitary(G, [1, 2], label=f'G{layer}_1')
    
    return qc
```
]

= Example Submission

#text(size: 9pt)[
```python
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate
from scipy.stats import unitary_group
import numpy as np

# Build target
U_target = build_mod_mult_unitary(a=7, N=15, n_qubits=4)

# Example: Brick-wall ansatz with optimized SU(4) gates
qc = QuantumCircuit(4)

# Layer 1: gates on (0,1) and (2,3)
G1 = unitary_group.rvs(4)  # Replace with optimized parameters
G2 = unitary_group.rvs(4)
qc.unitary(G1, [0, 1], label='G1')
qc.unitary(G2, [2, 3], label='G2')

# Layer 2: gate on (1,2)
G3 = unitary_group.rvs(4)
qc.unitary(G3, [1, 2], label='G3')

# Layer 3: gates on (0,1) and (2,3)
G4 = unitary_group.rvs(4)
G5 = unitary_group.rvs(4)
qc.unitary(G4, [0, 1], label='G4')
qc.unitary(G5, [2, 3], label='G5')

# Verify
result = verify_circuit(qc, U_target)
print(f"Fidelity: {result['fidelity']:.6f}")
print(f"SU(4) count: {result['su4_count']}")
print(f"Valid (F >= 0.99): {result['is_valid']}")

# Visualize circuit
print(qc.draw())
```
]

= Qiskit Circuit Visualization

The submitted circuit can be visualized using Qiskit:

#text(size: 9pt)[
```python
# Draw circuit
from qiskit.visualization import circuit_drawer
circuit_drawer(qc, output='mpl', filename='circuit.png')

# Or in text format:
# q_0: ─[G1]─────[G4]─
#      │        │
# q_1: ─┘─[G3]──┘─────
#         │
# q_2: ───┴───[G2][G5]
#             │   │
# q_3: ───────┘───┘───
```
]

= Baselines

#table(
  columns: (1fr, auto, auto),
  inset: 6pt,
  align: (left, center, left),
  [*Method*], [*SU(4)*], [*Notes*],
  [Generic 4-qubit permutation], [$<= 6$], [Upper bound #cite(<shende2004>)],
  [Random initialization], [6], [Brick-wall, 2 layers],
  [Theoretical lower bound], [$>= 3$], [Connectivity constraint],
)

*Challenge*: Can you achieve $F >= 0.99$ with fewer than 6 SU(4) gates?

= Evaluation

+ *Primary*: SU(4) gate count
+ *Secondary*: Circuit depth, achieved fidelity

#v(1em)
#line(length: 100%, stroke: 0.3pt)
#v(0.5em)

#bibliography("refs.bib", title: "References", style: "ieee")
