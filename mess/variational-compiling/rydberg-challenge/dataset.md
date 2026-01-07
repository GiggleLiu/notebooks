# Dataset Description

**File**: `rydberg_challenge.npz`

150 target unitaries of varying sizes. **The optimal pulse count is unknown** for most targets.

## Target Categories

| ID | Qubits | Type | Description |
|----|--------|------|-------------|
| 1-20 | 2 | Clifford | CZ, CNOT, CPhase variants |
| 21-40 | 2 | Parametric | Two-RZZ constructions |
| 41-70 | 2 | Random | Haar-random SU(4) |
| 71-85 | 2 | SWAP-like | SWAP, √SWAP, etc. |
| 86-100 | 2 | Physical | Ising evolution |
| 101-110 | 3 | Trotterized | Heisenberg/Ising evolution |
| 111-120 | 3 | Random | Haar-random SU(8) |
| 121-130 | 4 | Trotterized | TFIM, XXZ models |
| 131-140 | 4 | Random | Haar-random SU(16) |
| 141-145 | 5 | Trotterized | Physical Hamiltonians |
| 146-150 | 6 | Trotterized | Large-scale evolution |

## Loading the Dataset

```python
import numpy as np

data = np.load('rydberg_challenge.npz', allow_pickle=True)
unitaries = data['unitaries']      # Array of target unitary matrices
n_qubits_list = data['n_qubits']   # Number of qubits for each target
types = data['types']              # Category label for each target

# Example: access target 0
U = unitaries[0]
n = n_qubits_list[0]
print(f"Target 1: {n}-qubit {types[0]} unitary, shape {U.shape}")
```

## Verification Program

```python
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, process_fidelity
from qiskit.converters import circuit_to_dag

def count_pulses(qc: QuantumCircuit) -> int:
    """
    Count pulses assuming same-type gates can be parallel.
    Each pulse contains either all 1-qubit or all 2-qubit gates.
    """
    dag = circuit_to_dag(qc)
    pulses = 0
    
    for layer in dag.layers():
        ops = list(layer['graph'].op_nodes())
        if not ops:
            continue
        
        # Check if layer has mixed gate types
        has_1q = any(node.op.num_qubits == 1 for node in ops)
        has_2q = any(node.op.num_qubits == 2 for node in ops)
        
        if has_1q and has_2q:
            pulses += 2  # Need separate pulses for 1q and 2q
        else:
            pulses += 1
    
    return pulses

def verify_circuit(qc: QuantumCircuit, U_target: np.ndarray) -> dict:
    """
    Verify a Qiskit circuit against the target unitary.
    
    Returns:
        dict with fidelity, pulse_count, rzz_count, is_valid
    """
    circuit_op = Operator(qc)
    target_op = Operator(U_target)
    fidelity = process_fidelity(circuit_op, target_op)
    
    rzz_count = sum(1 for inst in qc.data if inst.operation.name == 'rzz')
    pulse_count = count_pulses(qc)
    
    return {
        'fidelity': fidelity,
        'infidelity': 1 - fidelity,
        'pulse_count': pulse_count,
        'rzz_count': rzz_count,
        'is_valid': (1 - fidelity) < 1e-8
    }
```

## Example Circuit

```python
def example_circuit():
    qc = QuantumCircuit(4)
    
    # Pulse 1: parallel single-qubit gates
    qc.ry(0.5, 0)
    qc.ry(0.3, 1)
    qc.ry(0.7, 2)
    qc.ry(0.2, 3)
    
    # Pulse 2: parallel RZZ gates on (0,1) and (2,3)
    qc.rzz(0.4, 0, 1)
    qc.rzz(0.5, 2, 3)
    
    # Pulse 3: single-qubit gates
    qc.rx(0.2, 1)
    qc.rx(0.1, 2)
    
    # Pulse 4: RZZ on (1,2)
    qc.rzz(0.6, 1, 2)
    
    return qc  # Total: 4 pulses
```
