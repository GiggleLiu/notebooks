# 里德堡原子门编译挑战

## 背景

里德堡原子阵列是一种领先的量子计算平台，其原生的两量子比特相互作用是 ZZ 门。本挑战的重点是以最少的脉冲数编译任意酉算符。

## 门集

- **单量子比特门**：任意 U(2) 旋转
- **双量子比特门**：$R_{ZZ}(\theta) = \exp(-i \frac{\theta}{2} Z \otimes Z)$
- **连通性**：任意（全连接）
- **并行化**：作用在不相交量子比特上的同类型门可在一个脉冲内执行

---

## 问题描述

> **给定**：一个酉算符 $U$  
> **求解**：具有最少脉冲数的电路分解

### 约束条件

- **允许的门**：任意量子比特对上的 RZZ(θ)，任意单量子比特门
- **精度要求**：不保真度 $< 10^{-8}$
- **并行性**：每个脉冲只能包含同类型的门（全部单量子比特门或全部双量子比特门）

### 目标

$$\min \text{脉冲数} \quad \text{满足} \quad 1 - F(U, V) < 10^{-8}$$

其中保真度 $F = |\text{Tr}(U^\dagger V)|^2 / d^2$，$d = 2^n$。

---

## 提交格式

```python
from qiskit import QuantumCircuit
import numpy as np

def compile(U_target: np.ndarray, n_qubits: int) -> QuantumCircuit:
    """
    将目标酉算符编译为 RZZ + 单量子比特门。
    
    参数:
        U_target: 目标酉矩阵，numpy 数组 (2^n × 2^n)
        n_qubits: 量子比特数
    
    返回:
        qc: 仅使用以下门的 QuantumCircuit：
            - 单量子比特门 (rx, ry, rz, u 等)
            - 任意量子比特对上的 rzz 门
    
    约束:
        - 目标精度：不保真度 < 1e-8
        - 最小化脉冲数（同类型并行门 = 1 个脉冲）
    """
    qc = QuantumCircuit(n_qubits)
    # 在此实现你的算法
    return qc
```
