# 数据集说明

**文件**：`rydberg_challenge.npz`

包含 150 个不同规模的目标酉算符。**大多数目标的最优脉冲数未知**。

## 目标分类

| 编号 | 量子比特数 | 类型 | 描述 |
|----|--------|------|-------------|
| 1-20 | 2 | Clifford 门 | CZ, CNOT, CPhase 变体 |
| 21-40 | 2 | 参数化门 | 两个 RZZ 构造 |
| 41-70 | 2 | 随机 | Haar 随机 SU(4) |
| 71-85 | 2 | SWAP 类 | SWAP, √SWAP 等 |
| 86-100 | 2 | 物理系统 | Ising 演化 |
| 101-110 | 3 | Trotter 化 | Heisenberg/Ising 演化 |
| 111-120 | 3 | 随机 | Haar 随机 SU(8) |
| 121-130 | 4 | Trotter 化 | TFIM, XXZ 模型 |
| 131-140 | 4 | 随机 | Haar 随机 SU(16) |
| 141-145 | 5 | Trotter 化 | 物理哈密顿量 |
| 146-150 | 6 | Trotter 化 | 大规模演化 |

## 加载数据集

```python
import numpy as np

data = np.load('rydberg_challenge.npz', allow_pickle=True)
unitaries = data['unitaries']      # 目标酉矩阵数组
n_qubits_list = data['n_qubits']   # 每个目标的量子比特数
types = data['types']              # 每个目标的类别标签

# 示例：访问目标 0
U = unitaries[0]
n = n_qubits_list[0]
print(f"目标 1: {n} 量子比特 {types[0]} 酉算符，形状 {U.shape}")
```

## 验证程序

```python
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, process_fidelity
from qiskit.converters import circuit_to_dag

def count_pulses(qc: QuantumCircuit) -> int:
    """
    计算脉冲数，假设同类型门可以并行。
    每个脉冲包含全部单量子比特门或全部双量子比特门。
    """
    dag = circuit_to_dag(qc)
    pulses = 0
    
    for layer in dag.layers():
        ops = list(layer['graph'].op_nodes())
        if not ops:
            continue
        
        # 检查该层是否有混合门类型
        has_1q = any(node.op.num_qubits == 1 for node in ops)
        has_2q = any(node.op.num_qubits == 2 for node in ops)
        
        if has_1q and has_2q:
            pulses += 2  # 需要分开的脉冲用于单量子比特门和双量子比特门
        else:
            pulses += 1
    
    return pulses

def verify_circuit(qc: QuantumCircuit, U_target: np.ndarray) -> dict:
    """
    验证 Qiskit 电路与目标酉算符的匹配程度。
    
    返回:
        包含 fidelity, pulse_count, rzz_count, is_valid 的字典
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

## 示例电路

```python
def example_circuit():
    qc = QuantumCircuit(4)
    
    # 脉冲 1: 并行单量子比特门
    qc.ry(0.5, 0)
    qc.ry(0.3, 1)
    qc.ry(0.7, 2)
    qc.ry(0.2, 3)
    
    # 脉冲 2: (0,1) 和 (2,3) 上的并行 RZZ 门
    qc.rzz(0.4, 0, 1)
    qc.rzz(0.5, 2, 3)
    
    # 脉冲 3: 单量子比特门
    qc.rx(0.2, 1)
    qc.rx(0.1, 2)
    
    # 脉冲 4: (1,2) 上的 RZZ 门
    qc.rzz(0.6, 1, 2)
    
    return qc  # 总计: 4 个脉冲
```
