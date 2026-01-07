# Variational Quantum Compilation Challenge

## Background

Variational quantum compilation (VQC) approximates a target unitary $U_{\text{target}}$ with a parameterized circuit $V(\boldsymbol{\theta})$ by maximizing the fidelity [1]:

$$F = \frac{|\text{Tr}(U_{\text{target}}^\dagger V(\boldsymbol{\theta}))|^2}{d^2}$$

---

## Shor's Algorithm: The Bottleneck

| Component | Complexity | Notes |
|-----------|------------|-------|
| Hadamard layer | $O(n)$ | Initial superposition |
| **Modular exponentiation** | **$O(n^3)$** | **Bottleneck: $2n$ controlled mod-mult** |
| Inverse QFT | $O(n^2)$ | $n(n-1)/2$ controlled rotations |

---

## Challenge

### Target Unitary

Compile the **4-qubit modular multiplication** from Shor's algorithm for factoring $N = 15$:

$$U_7: |y\rangle \rightarrow |7 \cdot y \mod 15\rangle$$

This is a **permutation matrix** with period $r = 4$ (since $7^4 \equiv 1 \mod 15$).

### Circuit Structure

```
|q0⟩ ──[G1]──────[G4]── |q'0⟩
      │          │
|q1⟩ ──┴──[G2][G3]┴──── |q'1⟩
           │  │
|q2⟩ ──────┴──┴──────── |q'2⟩

|q3⟩ ────────────────── |q'3⟩
```

Each $G_i \in \text{SU}(4)$ acts on adjacent qubits. Connectivity: linear chain.

### Objective

> **Goal:** Find circuit $V(\boldsymbol{\theta})$ with **minimum SU(4) count** achieving $F \geq 0.99$

---

## Target Specification

### Permutation Table

| $y$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|-----|---|---|---|---|---|---|---|---|
| $7y \mod 15$ | 0 | 7 | 14 | 6 | 13 | 5 | 12 | 4 |

| $y$ | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|-----|---|---|---|---|---|---|---|---|
| $7y \mod 15$ | 11 | 3 | 10 | 2 | 9 | 1 | 8 | 15 |

**Key cycle:** $1 \rightarrow 7 \rightarrow 4 \rightarrow 13 \rightarrow 1$ (period 4)

### Target Unitary (Qiskit)

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

---

## Verification Program (Qiskit)

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
    
    # Count custom unitaries (SU4 gates)
    su4_count = sum(1 for inst in qc.data 
                    if isinstance(inst.operation, UnitaryGate) 
                    and inst.operation.num_qubits == 2)
    
    return {
        "fidelity": fidelity,
        "su4_count": su4_count,
        "is_valid": fidelity >= 0.99
    }

def build_ansatz(n_layers: int = 2) -> QuantumCircuit:
    """Build a brick-wall ansatz with SU(4) gates."""
    qc = QuantumCircuit(4)
    
    for layer in range(n_layers):
        if layer % 2 == 0:
            for q in [0, 2]:
                if q + 1 < 4:
                    from scipy.stats import unitary_group
                    G = unitary_group.rvs(4)
                    qc.unitary(G, [q, q+1], label=f'G{layer}_{q}')
        else:
            from scipy.stats import unitary_group
            G = unitary_group.rvs(4)
            qc.unitary(G, [1, 2], label=f'G{layer}_1')
    
    return qc
```

---

## Example Submission

```python
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate
from scipy.stats import unitary_group
import numpy as np

# Build target
U_target = build_mod_mult_unitary(a=7, N=15, n_qubits=4)

# Example: Brick-wall ansatz with SU(4) gates
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

---

## Baselines

| Method | SU(4) Count | Notes |
|--------|-------------|-------|
| Generic 4-qubit permutation | ≤6 | Upper bound [6] |
| Random initialization | 6 | Brick-wall, 2 layers |
| Theoretical lower bound | ≥3 | Connectivity constraint |

**Challenge:** Can you achieve F ≥ 0.99 with fewer than 6 SU(4) gates?

---

## Evaluation

1. **Primary:** SU(4) gate count
2. **Secondary:** Circuit depth, achieved fidelity

---

## References

[1] S. Khatri et al., "Quantum-assisted quantum compiling," *Quantum* **3**, 140 (2019).

[2] A. Botea et al., "On the complexity of quantum circuit compilation," *Int. Symp. Combinatorial Search* (2018).

[3] S. Beauregard, "Circuit for Shor's algorithm using 2n+3 qubits," *Quantum Inf. Comput.* **3**(2), 175-185 (2003).

[4] L. M. K. Vandersypen et al., "Experimental realization of Shor's quantum factoring algorithm using nuclear magnetic resonance," *Nature* **414**, 883-887 (2001).

[5] N. Khaneja and S. J. Glaser, "Cartan decomposition of SU(2^n)," *Chem. Phys.* **267**, 11-23 (2001).

[6] V. V. Shende et al., "Minimal universal two-qubit controlled-NOT-based circuits," *Phys. Rev. A* **69**, 062321 (2004).
