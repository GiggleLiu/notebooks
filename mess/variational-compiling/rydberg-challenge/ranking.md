# Ranking & Evaluation

## Scoring

| Metric | Description |
|--------|-------------|
| **Primary** | Total pulse count (lower is better) |
| **Constraint** | All infidelities must be < 10⁻⁸ |

Submissions that fail the precision constraint on any target are **disqualified**.

---

## Evaluation Script

```python
import numpy as np
from qiskit.quantum_info import Operator

# Load dataset
data = np.load('rydberg_challenge.npz', allow_pickle=True)
unitaries = data['unitaries']
n_qubits_list = data['n_qubits']
types = data['types']

# Evaluate submission
total_pulses = 0
failed = 0
results = []

for i, (U, n) in enumerate(zip(unitaries, n_qubits_list)):
    qc = compile(U, n)  # Your compile function
    result = verify_circuit(qc, U)
    
    total_pulses += result['pulse_count']
    results.append(result)
    
    if not result['is_valid']:
        failed += 1
        print(f"FAIL Target {i+1}: infidelity={result['infidelity']:.2e}")

print(f"\n{'='*40}")
print(f"Total pulses: {total_pulses}")
print(f"Failed targets: {failed}/150")
print(f"{'='*40}")

if failed == 0:
    print(f"\n✓ SCORE: {total_pulses} pulses")
else:
    print(f"\n✗ DISQUALIFIED: {failed} targets failed precision requirement")
```

---

## Baseline Estimates

| Qubits | Targets | Naïve Pulses | KAK-Based |
|--------|---------|--------------|-----------|
| 2 | 100 | ~1400 | ~700 |
| 3 | 20 | ~600 | ~300 |
| 4 | 20 | ~1200 | ~400 |
| 5 | 5 | ~500 | ~150 |
| 6 | 5 | ~800 | ~200 |
| **Total** | **150** | **~4500** | **~1750** |

### Baseline Methods

- **Naïve**: Standard Trotterization with sequential gates
- **KAK-Based**: Uses Cartan decomposition for 2-qubit blocks, then optimal scheduling

---

## Leaderboard

| Rank | Team | Total Pulses | Notes |
|------|------|--------------|-------|
| - | Baseline (Naïve) | ~4500 | Sequential execution |
| - | Baseline (KAK) | ~1750 | Optimal 2-qubit decomposition |
| 1 | *Your team here* | ? | |

---

*Can you minimize the total pulse count with efficient gate scheduling?*
