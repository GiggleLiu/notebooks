#set text(font: "New Computer Modern", size: 11pt)
#set page(margin: 1in)
#set par(justify: true)
#set heading(numbering: "1.1")

// Bra-ket notation
#let ket(x) = $lr(|#x chevron.r)$
#let bra(x) = $lr(chevron.l #x|)$
#let braket(x, y) = $lr(chevron.l #x|#y chevron.r)$

#align(center)[
  #text(size: 18pt, weight: "bold")[Tensor Network Methods for Variational Quantum Compilation]
  
  #v(0.5em)
  #text(size: 11pt, style: "italic")[Literature Survey -- January 2026]
]

#v(1em)

= Introduction

Variational quantum compilation (VQC) seeks to approximate target quantum operations with shallow, hardware-executable circuits through optimization. The integration of tensor network methods has emerged as one of the most powerful approaches in this field, providing:
- Efficient classical simulation and cost evaluation
- Physics-informed ansatz design
- Hybrid classical-quantum optimization strategies
- Theoretical tools for understanding entanglement and approximation errors

This survey focuses on the synergy between tensor networks and variational quantum compilation.

= Background: Variational Quantum Compilation

== Problem Formulation

Given a target unitary $U_"target"$, VQC finds parameters $bold(theta)$ for a parameterized ansatz $V(bold(theta))$ that minimizes the infidelity:
$ C = 1 - |"Tr"(U_"target"^dagger V(bold(theta)))|^2 / d^2 $
where $d = 2^n$ is the Hilbert space dimension.

== Key Methods

- *QAQC* #cite(<khatri2019qaqc>, form: "normal"): Foundational framework using quantum hardware for cost evaluation
- *LVQC* #cite(<mizuta2022local>, form: "normal"): Local cost functions via Lieb-Robinson bounds for scalability
- *VFF* #cite(<commeau2020vff>, form: "normal"): Variational fast forwarding for time evolution compression

= Tensor Networks as Circuit Representations

== Circuit-Tensor Network Correspondence

Quantum circuits naturally admit tensor network representations. Each gate is a tensor, and qubit lines are contracted indices:
$ U = sum_(i_1, ..., i_k) G_1^(i_1) G_2^(i_2) dots.c G_k^(i_k) $

This correspondence is foundational to tensor network methods in VQC, enabling:
- *Classical simulation*: Contract tensor networks for circuits with bounded entanglement #cite(<vidal2003efficient>, form: "normal")
- *Fidelity computation*: Evaluate $|"Tr"(U^dagger V)|^2$ via contraction #cite(<markov2008simulating>, form: "normal")
- *Gradient estimation*: Compute parameter gradients classically

== Contraction Complexity

The cost of tensor network contraction depends on:
- *Treewidth*: Exponential in the treewidth of the circuit graph
- *Bond dimension*: Intermediate tensor sizes during contraction
- *Contraction ordering*: Optimal ordering is NP-hard; heuristics are essential #cite(<gray2021hyper>, form: "normal")

For circuits with treewidth $w$, contraction scales as $O(2^w)$, making tensor methods efficient for structured circuits.

= Tensor Network Ansätze for VQC

== Matrix Product States (MPS)

MPS-inspired circuits use sequential two-qubit gates in brick-wall patterns:
$ #ket($psi$) = sum_(s_1, ..., s_n) A_1^(s_1) A_2^(s_2) dots.c A_n^(s_n) #ket($s_1 dots.c s_n$) $

Properties:
- Bond dimension $chi$ controls entanglement: $S_"max" = log chi$
- Area-law states (1D ground states) are efficiently representable
- Natural mapping to linear qubit connectivity #cite(<schollwock2011dmrg>, form: "normal")

== Tree Tensor Networks (TTN)

TTN employ hierarchical tree structures:
- Logarithmic depth: $O(log n)$ layers for $n$ qubits
- Capture hierarchical correlations efficiently
- Well-suited for systems with multi-scale structure

== MERA Circuits

Multi-scale Entanglement Renormalization Ansatz #cite(<evenbly2009algorithms>, form: "normal") interleaves:
- *Disentanglers*: Remove short-range entanglement
- *Isometries*: Coarse-grain degrees of freedom

MERA can represent critical systems with logarithmic entanglement scaling, beyond MPS capabilities.

== Quantum Circuit Tensor Networks (QCTN)

#cite(<qctn2020>, form: "prose") introduced QCTN, replacing dense unitaries with parameterized local circuits:
$ T_"QCTN" = product_(l=1)^L U_l (bold(theta)_l) $

Advantages:
- Exponentially more expressive than classical tensor networks at fixed parameter count
- Avoid $O(chi^4)$ scaling of dense unitaries
- Trainable via quantum gradient methods

Demonstrated applications:
- Many-body ground state preparation with high fidelity
- Efficient representation beyond area-law states

= Approximate Quantum Compilation with Tensor Networks

== AQC-Tensor Framework

IBM's AQC-Tensor #cite(<aqc_tensor>, form: "normal") provides a systematic workflow:

+ *Target representation*: Express $U_"target"$ as a tensor network (e.g., from Trotter expansion)
+ *Ansatz construction*: Build parameterized circuit matching target structure
+ *Cost evaluation*: Use tensor contraction for efficient fidelity computation
+ *Optimization*: Gradient descent with classically computed gradients
+ *Compression*: Iteratively reduce depth while maintaining fidelity

Particularly effective for:
- Time evolution: Compress deep Trotter circuits to shallow approximations
- State preparation: Reduce depth for structured initial states
- Hamiltonian simulation: Approximate $e^(-i H t)$ within hardware constraints

== Tensor Network Contraction for Cost Functions

For moderate-entanglement targets, the Hilbert-Schmidt cost:
$ C = 1 - |"Tr"(U_"target"^dagger V(bold(theta)))|^2 / d^2 $
can be computed via tensor contraction without quantum hardware.

This enables:
- *Pre-screening*: Test ansatz expressibility before quantum execution
- *Gradient computation*: Full gradient via automatic differentiation
- *Hyperparameter tuning*: Optimize circuit structure classically

== Iterative Circuit Compression

Layer-by-layer compression algorithm:
+ Represent first $k$ circuit layers as tensor network
+ Variationally optimize $(k-1)$-layer approximation
+ Iterate until target depth achieved
+ Fine-tune on quantum hardware if needed

Error accumulation can be bounded via singular value truncation analysis.

= Hybrid Classical-Quantum Strategies

== Classical Pre-optimization

Tensor networks enable efficient warm-starting:
+ Simulate target circuit as tensor network
+ Optimize ansatz parameters classically
+ Transfer parameters to quantum device
+ Refine where classical simulation fails

This dramatically reduces quantum resource requirements for compilation.

== Entanglement-Guided Ansatz Design

Tensor network analysis reveals entanglement structure:
- Compute entanglement entropy across cuts
- Identify bond dimensions needed for target fidelity
- Design ansatz with matching entanglement capacity

For area-law targets, MPS-inspired circuits suffice; volume-law requires deeper circuits or QCTN.

== Tensor Network Feature Extraction

#cite(<chen2020hybrid>, form: "prose") demonstrated hybrid architectures:
- MPS compresses classical input data
- Variational quantum circuit processes compressed features
- End-to-end training via backpropagation

This leverages tensor networks for dimensionality reduction before quantum processing.

= Software Frameworks

== TensorCircuit

#cite(<tensorcircuit2022>, form: "prose") provides a comprehensive framework:
- Tensor network-based circuit simulation
- Automatic differentiation (JAX, TensorFlow, PyTorch backends)
- JIT compilation and GPU/TPU acceleration
- Noise simulation capabilities
- Monte Carlo and MPS simulation modes

Key for large-scale VQC research with classical optimization.

== Additional Tools

- *Quimb*: Python tensor network library with circuit simulation support
- *ITensor*: High-performance tensor networks (Julia/C++)
- *TeNPy*: Tensor networks for many-body physics
- *cotengra*: Optimal contraction path finding #cite(<gray2021hyper>, form: "normal")

= Theoretical Insights from Tensor Networks

== Expressibility Analysis

Tensor networks quantify circuit expressibility:
- MPS circuits can represent states with entanglement $S <= log chi$
- MERA captures logarithmic corrections
- Random circuits approach Haar-random at sufficient depth

This guides ansatz selection for specific compilation targets.

== Approximation Bounds

Singular value decomposition provides error bounds:
$ ||U - tilde(U)||_F <= sum_i sigma_i^"trunc" $

Truncating small singular values during compression yields controlled approximation errors.

== Barren Plateau Mitigation

Tensor network structure can mitigate barren plateaus:
- Local cost functions from tensor network decomposition
- Layer-wise training inspired by DMRG sweeps
- Entanglement-limited ansätze avoid exponential gradient decay

= Minimal Tensor Network Representations

A fundamental question in tensor network methods is: _given a quantum operator, what is its minimal tensor network representation?_ This question has multiple interpretations and remains partially open.

== Notions of Minimality

"Minimal" can refer to several measures:

+ *Bond dimension*: Minimize max bond dimension $chi_"max"$ across all edges
+ *Total parameters*: Minimize $sum_e chi_e$ or total tensor entries
+ *Network structure*: Find optimal graph topology (MPS, TTN, PEPS, etc.)
+ *Approximation trade-off*: Minimize bond dimension for fixed error $epsilon$

These objectives often conflict---optimal structure may differ from optimal bond dimension within a fixed structure.

== Approaches for Finding Minimal Representations

=== SVD-Based Sequential Decomposition

The most direct approach uses successive singular value decompositions:

+ Reshape operator as matrix: $O_(i_1 ... i_k, j_1 ... j_l)$
+ Apply SVD: $O = U Sigma V^dagger$
+ Truncate to rank $chi$: keep largest $chi$ singular values
+ Reshape $U$ as first tensor, continue with $Sigma V^dagger$

For MPS/MPO, this yields the *TT-SVD* (tensor-train SVD) algorithm with complexity $O(n d^2 chi^2)$.

*Guarantees*: TT-SVD is quasi-optimal---the resulting bond dimension is at most a factor of $sqrt(n)$ larger than optimal for the same approximation error.

=== Variational Optimization (DMRG-style)

Iteratively optimize each tensor while keeping others fixed:

+ Initialize with random or SVD-based guess
+ Sweep through tensors, optimizing each via local eigenvalue problem
+ Repeat until convergence

*Advantages*: Can achieve lower bond dimensions than SVD for structured operators
*Disadvantages*: No global optimality guarantee; may converge to local minima

=== Density Matrix Approach

For operators arising from physical systems:
+ Compute reduced density matrices for subsystems
+ Bond dimension bounded by rank of reduced density matrix
+ Entanglement entropy provides lower bound: $chi >= 2^(S)$

=== Automatic MPO Construction

For Hamiltonians as sums of local terms $H = sum_j h_j$:
+ Finite-state automaton construction
+ Build MPO tensors encoding operator strings
+ Bond dimension scales with number of distinct operator types

This gives *exact* (not approximate) MPO with minimal bond dimension for the given structure.

=== Rank-Revealing Decompositions

Beyond SVD, other decompositions can find minimal representations:
- *Interpolative decomposition*: Select actual columns/rows as basis
- *Randomized methods*: Efficient for large-scale problems
- *Cross approximation*: Adaptive sampling for low-rank structure

== Complexity of Finding Optimal Representations

=== Tensor Rank is NP-Hard

Computing the exact tensor rank (minimal CP decomposition) is NP-hard. Related problems:

- *Border rank*: Limit of rank-$r$ approximations (also hard)
- *Multilinear rank*: Tuple of matrix ranks (polynomial time via SVD)
- *Tree tensor network rank*: Depends on tree structure

=== Structure Optimization

Finding the optimal tensor network _structure_ (graph topology) is generally intractable:

- Exponentially many possible graphs for $n$ tensors
- Treewidth optimization is NP-complete
- Heuristics: greedy algorithms, simulated annealing, genetic algorithms

=== Approximation Hardness

Even approximate rank computation is hard:
- No PTAS (polynomial-time approximation scheme) unless P=NP
- Best known approximation: $O(n^(1-epsilon))$ factor for rank-$r$ tensors

== Open Questions

=== Optimal Bond Dimension Certification

*Problem*: Given operator $O$ and target fidelity $1-epsilon$, what is the minimal bond dimension $chi^*$?

- Lower bounds from entanglement entropy are often loose
- No efficient algorithm to certify optimality
- Gap between achievable and provable bounds

=== Global vs. Local Minima in Variational Methods

*Problem*: When does variational optimization (DMRG) find the global minimum?

- Known to work well for gapped 1D systems
- Can fail for critical systems or 2D
- Landscape analysis largely unexplored

=== Extension to 2D and Beyond

*Problem*: Generalize minimal representation theory to PEPS and higher-dimensional networks.

- PEPS contraction is \#P-hard
- No analog of Schmidt decomposition in 2D
- Approximate methods (boundary MPS) introduce uncontrolled errors

=== Long-Range Interactions

*Problem*: Efficient representation of operators with power-law decay $1\/r^alpha$.

- MPO bond dimension can grow polynomially with system size
- Exponential representations exist but may not be minimal
- Trade-off between locality and bond dimension

=== Operator Entanglement vs. State Entanglement

*Problem*: Relate operator entanglement (MPO bond dimension) to state entanglement under evolution.

- Time evolution can increase operator entanglement linearly
- Implications for Hamiltonian simulation complexity
- Connection to quantum chaos and scrambling

=== Certifiable Approximation Bounds

*Problem*: Given approximate TN representation, bound the approximation error efficiently.

- SVD truncation gives Frobenius norm bound
- Spectral norm bounds harder to obtain
- Application-specific error metrics (e.g., expectation values)

== Practical Recommendations

For practitioners seeking minimal representations:

+ *Start with TT-SVD*: Provides good baseline with guarantees
+ *Refine with variational*: DMRG sweeps can reduce bond dimension
+ *Exploit structure*: Use automatic MPO for sum-of-local-terms
+ *Monitor convergence*: Track singular value spectrum for truncation decisions
+ *Consider application*: Minimal for storage may differ from minimal for contraction

= Applications

== Hamiltonian Simulation

Primary application of tensor network VQC:
- Compress Trotter sequences: $e^(-i H t) approx (e^(-i H_1 Delta t) e^(-i H_2 Delta t))^(t\/Delta t)$
- Reduce from $O(t\/epsilon)$ to $O(1)$ depth via variational compression
- Tensor network contraction evaluates fidelity efficiently

== Quantum Chemistry

Compilation of molecular simulation circuits:
- UCC ansätze have tensor network structure
- Compress excitation operators to hardware-native gates
- Demonstrated 75%+ CNOT reduction in benchmarks #cite(<cowtan2020compilation>, form: "normal")

== State Preparation

Tensor network methods for initial state compilation:
- Ground states of local Hamiltonians have area-law entanglement
- MPS-inspired circuits achieve efficient preparation
- QCTN extends to more entangled states

= Challenges and Open Problems

== Entanglement Barriers

Highly entangled targets require:
- Large bond dimensions (exponential classical cost)
- Deep circuits (hardware limitations)
- Hybrid approaches may help but don't fully resolve

== Contraction Scaling

For general circuits:
- Treewidth grows with circuit depth and connectivity
- Exact contraction becomes intractable
- Approximate contraction methods introduce errors

== Hardware Topology Mismatch

Tensor network structures may not match device connectivity:
- MPS natural for linear chains
- 2D hardware may prefer PEPS-like structures
- SWAP overhead can negate compression benefits

= Future Directions

+ *Optimal contraction for VQC*: Develop contraction algorithms tailored to variational compilation workloads
+ *Entanglement-aware compilation*: Automatically design ansätze matching target entanglement profiles
+ *Error-bounded compression*: Provide rigorous fidelity guarantees for compiled circuits
+ *Hardware-topology-aware tensor networks*: Design tensor network structures matching specific quantum processors
+ *Scalable hybrid methods*: Extend tensor network VQC to 100+ qubit systems

= Conclusion

Tensor network methods have become indispensable for variational quantum compilation, providing:
- Efficient classical simulation for cost evaluation and gradient computation
- Physics-informed ansatz design through MPS, TTN, MERA, and QCTN structures
- Hybrid optimization strategies that minimize quantum resource requirements
- Theoretical tools for understanding expressibility and approximation errors

Software frameworks like TensorCircuit and AQC-Tensor make these methods accessible for practical VQC research. As quantum hardware scales, the synergy between tensor networks and variational compilation will be essential for bridging the gap between theoretical algorithms and executable circuits.

Key challenges remain in handling highly entangled targets and matching tensor network structures to hardware topologies, but ongoing research continues to advance this productive intersection of quantum computing and tensor network theory.

#pagebreak()

#bibliography("refs.bib", title: "References", style: "ieee")
