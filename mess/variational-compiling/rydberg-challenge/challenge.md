# Rydberg Atom Gate Compilation Challenge

## Background

Rydberg atom arrays are a leading quantum computing platform where the native two-qubit interaction is the ZZ gate. This challenge focuses on compiling arbitrary unitaries with minimal number of pulses.

## Gate Set

- **Single-qubit**: Arbitrary U(2) rotation
- **Two-qubit**: $R_{ZZ}(\theta) = \exp(-i \frac{\theta}{2} Z \otimes Z)$
- **Connectivity**: Arbitrary (all-to-all)
- **Parallelization**: Same-type gates on disjoint qubits execute in one pulse

---

## Problem Statement

> **Given**: A unitary $U$  
> **Find**: Circuit decomposition with minimum number of pulses

### Constraints

- **Allowed gates**: RZZ(θ) on any pair, arbitrary single-qubit gates
- **Precision**: Infidelity $< 10^{-8}$
- **Parallelism**: Only same-type gates can be parallel (all-1q or all-2q per pulse)

### Objective

$$\min \text{Pulse Count} \quad \text{such that} \quad 1 - F(U, V) < 10^{-8}$$

where fidelity $F = |\text{Tr}(U^\dagger V)|^2 / d^2$ and $d = 2^n$.

---

## Submission Format

```python
from qiskit import QuantumCircuit
import numpy as np

def compile(U_target: np.ndarray, n_qubits: int) -> QuantumCircuit:
    """
    Compile a target unitary to RZZ + single-qubit gates.
    
    Args:
        U_target: Target unitary as numpy array (2^n × 2^n)
        n_qubits: Number of qubits
    
    Returns:
        qc: QuantumCircuit using only:
            - Single-qubit gates (rx, ry, rz, u, etc.)
            - rzz gates on ANY qubit pair
    
    Constraints:
        - Target precision: infidelity < 1e-8
        - Minimize pulse count (same-type parallel gates = 1 pulse)
    """
    qc = QuantumCircuit(n_qubits)
    # Your implementation here
    return qc
```
