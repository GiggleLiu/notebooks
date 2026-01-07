# 排名与评估

## 评分标准

| 指标 | 描述 |
|--------|-------------|
| **主要指标** | 总脉冲数（越低越好） |
| **约束条件** | 所有不保真度必须 < 10⁻⁸ |

任何目标未能满足精度约束的提交将被**取消资格**。

---

## 评估脚本

```python
import numpy as np
from qiskit.quantum_info import Operator

# 加载数据集
data = np.load('rydberg_challenge.npz', allow_pickle=True)
unitaries = data['unitaries']
n_qubits_list = data['n_qubits']
types = data['types']

# 评估提交
total_pulses = 0
failed = 0
results = []

for i, (U, n) in enumerate(zip(unitaries, n_qubits_list)):
    qc = compile(U, n)  # 你的编译函数
    result = verify_circuit(qc, U)
    
    total_pulses += result['pulse_count']
    results.append(result)
    
    if not result['is_valid']:
        failed += 1
        print(f"失败 目标 {i+1}: 不保真度={result['infidelity']:.2e}")

print(f"\n{'='*40}")
print(f"总脉冲数: {total_pulses}")
print(f"失败目标: {failed}/150")
print(f"{'='*40}")

if failed == 0:
    print(f"\n✓ 得分: {total_pulses} 个脉冲")
else:
    print(f"\n✗ 取消资格: {failed} 个目标未达到精度要求")
```

---

## 基准估计

| 量子比特数 | 目标数 | 朴素方法脉冲数 | 基于 KAK 方法 |
|--------|---------|--------------|-----------|
| 2 | 100 | ~1400 | ~700 |
| 3 | 20 | ~600 | ~300 |
| 4 | 20 | ~1200 | ~400 |
| 5 | 5 | ~500 | ~150 |
| 6 | 5 | ~800 | ~200 |
| **总计** | **150** | **~4500** | **~1750** |

### 基准方法

- **朴素方法**：使用顺序门的标准 Trotter 化
- **基于 KAK 方法**：对双量子比特块使用 Cartan 分解，然后进行最优调度

---

## 排行榜

| 排名 | 团队 | 总脉冲数 | 备注 |
|------|------|--------------|-------|
| - | 基准（朴素） | ~4500 | 顺序执行 |
| - | 基准（KAK） | ~1750 | 最优双量子比特分解 |
| 1 | *你的团队* | ? | |

---

*你能通过高效的门调度来最小化总脉冲数吗？*
