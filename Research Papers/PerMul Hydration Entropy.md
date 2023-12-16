### Previously Used Computational Methods:

1. **Thermodynamic Integration (TI):**
   - **Description:** TI is a method used to calculate solvation entropy by gradually transforming a system from one state to another while measuring changes in thermodynamic properties.
   - **Shortcomings:** It lacks spatial resolution, meaning it doesn't provide detailed information about where changes are occurring. Additionally, it requires extensive sampling, making it challenging for complex systems.

2. **3D-two-phase-thermodynamics (3D-2PT):**
   - **Description:** A voxel-based approach that provides spatially resolved hydration entropies by treating the system as a superposition of gas-like and solid-like components.
   - **Shortcomings:** It relies on the assumption that the system can be treated as a combination of two phases, which might not be accurate for heterogeneous systems.

3. **Grid Cell Theory (GCT):**
   - **Description:** GCT provides solvent entropies by relying on a generalized Pauling’s residual ice entropy model and represents the system using a grid-based approach.
   - **Shortcomings:** It assumes a specific model for the rotational entropy term, which may not be universally applicable, and the grid-based representation might not capture all the nuances of molecular interactions.

4. **Grid Inhomogeneous Solvation Theory (GIST):**
   - **Description:** GIST approximates entropy using a truncated expansion of single-body and multibody correlation functions calculated on a three-dimensional grid.
   - **Shortcomings:** The expansion is usually truncated, missing some multibody effects. It might not fully capture complex molecular interactions.

5. **Per|Mut Method (Permutation Reduction and Mutual Information Expansion):**
   - **Description:** Per|Mut utilizes permutation reduction to enhance sampling efficiency and a mutual information expansion for spatially resolved hydration entropies. It directly samples the configuration space probability density.
   - **Shortcomings:** It addresses the sampling problem but still relies on approximations in calculating mutual information expansion. The method is based on particle positions rather than voxels.

### Definitions:

1. Quaternion Space:
Quaternions are mathematical entities that extend complex numbers. They are often used to represent orientations in three-dimensional space, making them suitable for describing rotational degrees of freedom in molecular dynamics simulations.

2. Composite Space:
The composite space mentioned in the statement likely involves a combination of quaternion space and Euclidean space. This composite space is used to represent the combined features relevant to the system under study.

3. Gibbs Factor:
The Gibbs factor, denoted as $N!$, arises from statistical mechanics and refers to the number of ways a system of $N$ indistinguishable particles (molecules, in this case) can be rearranged without changing the macroscopic properties of the system. In the context of permutation reduction, it is used to account for the fact that physically identical microstates, differing only by the permutation of indistinguishable water molecules, are counted redundantly. The Gibbs factor is essential for reducing the configuration space volume, making the sampling more efficient.

### Formulaes:

The total entropy is split into translation, rotation and a mutual information term $I_{trans-rot}$.

1. **Translational Entropy:**

   $S_{\text{trans}} = -k_B \int \frac{\text{d}p^N \text{d}x^N}{h^{3\text{N}}} $
    - $k_B$ : Boltzmann Constant
    - $h$ : Planks Constant
    - $Q$ : Normalized and Dimensionless Space Density
    - $Q = Z^{-1} exp [- \frac{H}{k_B T}]$
    - $H$ : Hamiltonian 
    - $Z$ : Partition Function

2. **Configurational Entropy:**

    $S_{\text{conf}} \approx \sum_{i=1}^{\text{N}} S_1(i) ~-~ \sum_{(j,k)\in \text{pairs}} I_2(j,k) ~+~  \sum_{(l,m,n)\in \text{\text{triples}}} I_3(l,m,n)$

3. **Permutation Reduction Distance:**

   $x_{r(i)}(t) = \sum_{i=1}^{N} \left\| \mathbf{r}_i(t) - \mathbf{r}_{\pi(i)} \right\|^2$

4. **Entropy Estimation:**

   $S_{\text{est}}^{(p)} = -k \sum_{i=1}^{n_f} \frac{1}{\psi(x_{r(i)}/r_i^{(k)})} \left[ \log(x_{r(i)}/r_i^{(k)}) - \psi(x_{r(i)}/r_i^{(k)}) \right]$
   where $p$ denotes the term (single, pair, or triple), $n_f$ is the number of configurations, $k$ is a fixed positive integer, $\psi$ is the digamma function, $r_i^{(k)}$ is the distance to the $k$-th neighbor, and $x_{r(i)}$ is the distance between configurations.

5. **Correlation Term:**

   $I_{\text{trans-rot}} \approx \sum_{(j, \tilde{k}) \in \text{pairs}} {I_2}(j,\tilde{k})$

   where $j$ denotes translational degrees of freedom, $\tilde{k}$ denotes rotational degrees of freedom of molecule $k$.

6. **Composite Metric:**

   $d(\mathbf{x}_1, \mathbf{x}_2) = \sqrt{d_{\text{eucl}}^2(\mathbf{x}_1, \mathbf{x}_2) + \xi \cdot d_{\text{quat}}^2(\mathbf{q}_1, \mathbf{q}_2)}$
   where $\xi$ is a scaling factor, $d_{\text{eucl}}$ is the Euclidean metric, $d_{\text{quat}}$ is the quaternion metric, and $\mathbf{q}_1, \mathbf{q}_2$ are quaternions.

7. **Volume of the Ball in Composite Metric Space:**

   $V_r = \int_0^{2\pi} \int_0^{\pi} \int_0^{R} 128 \sin^2(\phi/2) \sin\theta \, d\phi \, d\theta \, dR$
   where $R$ is the radius of the ball, $\phi$ and $\theta$ are spherical coordinates.

8. Pairwise and Triplewise Correlation Formulas

    - $I_2(j,k) = S_1(j) + S_1(k) - S_2(j,k)$
    - $I_3(l,m,n) = S_1(l) + S_1(m) + S_1(n) - S_2(l,m) - S_2(m,n) - S_2(l,n) + S_3(l,m,n)$

