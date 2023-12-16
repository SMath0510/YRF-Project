## Introduction

- **Simulation Techniques:**
  - Molecular Dynamics (MD) and Kinetic Monte Carlo (KMC) simulations are used to study ion conduction in solid materials.
  - KMC enables the simulation of different time scale moves by selecting events based on probabilities.

- **Graph Representation:**
  - Ion conduction is represented as movement through a graph, where vertices correspond to energy structures and edges connect based on transition abilities.
  - The graph includes periodic boundary conditions for realistic large system simulations.

- **Centrality Measures Adaptation:**
  - Common network centrality measures (degree, closeness, betweenness) are adapted to focus on ion conduction time instead of the number of steps.

- **Transformation of Measures:**
  - Degree becomes adjacent flow-IN and adjacent flow-OUT, indicating rates of flow to and from a site.
  - Closeness centrality is transformed into return-flow centrality, measuring how easily ions flow to other vertices and return.
  - Betweenness centrality is renamed flow-through centrality, characterizing a sit

- **Robustness Analysis:**
  - The robustness of return-flow and flow-through centrality measures is tested by introducing random noise to rate constants.

## Methods and Definitions

### Centrality Measures Equations:

#### Definitions:
- $A_i$ is the adjacency list of $i$
- $m_{ij}$ is the mean time to go from $i$ to $j$.
- $m_{ijk}$ is the mean time to go from $i$ to $j$ and then from $j$ to $k$.
- $k_{nl}$ is the rate constant for motion from site $n$ to $l$.

1. **Adjacent Flow-IN and Flow-OUT:**
   - $C_i^{adj, IN} = \sum_{t \in A_i} k_{t \to i}$
   - $C_i^{adj, OUT} = \sum_{t \in A_i} k_{i \to t}$

2. **Return-Flow Centrality:**
   - $C_i^{return-flow} = \frac{1}{\sum_j \frac{1}{m_{ij} + m_{ji}}}$

3. **Flow-through Centrality (Periodic Paths):**
   - $C_i^{flow-through} = \sum_{t} \prod_{j=1}^{N} P_{i_j \to i_{j+1}}$
  
4. **Flow-through Centrality (All Paths):**
   - $C_i^{flow-through} = \frac{1}{\sum_{ij} \frac{m_{ij}}{m_{itj}}}$
   
5. **Formulaes:**
    - $m_{ij} = \frac{Z_{jj} - Z_{ij}}{\pi_j} \sum_{n} \pi_n c_n + \sum_n (Z_{in} - Z_{jn})c_n$
    - $Z = (I-P+W)^{-1}$. ($Z = $ fundamental matrix for ergodic chains)

### Additional Key Points:

- **Rate Constants ($k_{ij}$):**
   - Calculated using harmonic transition-state theory, dependent on energy differences and temperature.

- **Normalized Probabilities ($p_{ij}$):**
   - Probability of moving from site $i$ to $j$ is normalized based on the sum of all possible transitions.

- **Essence:**
   - **Periodic Long-Range Paths:** Focuses on probabilities of ion movement along specific paths, considering periodic boundary conditions.
   - **All Paths:** Considers mean time to go between sites for all possible paths, providing insights into how each site influences overall ion conduction.
   - The inverse of the average ratio in the formulas determines the centrality, highlighting sites that either accelerate or decelerate ion flow between two specified sites.

### Conclusion:

1. **Centrality Measures Adaptation:**
   - Three common centrality measures in network theory are modified to focus on time rather than the number of steps in the context of ion conduction.
   - These measures are applied to study proton conduction in yttrium-doped barium zirconate.

2. **Robustness and Vertex Removal:**
   - Flow-Through Centrality exhibits robustness to the removal of the most central vertices, providing consistent results even in scenarios where the highest centrality site is blocked.
   - When considering all possible paths, Flow-Through Centrality highlights regions for both traps and highways, offering insights into probable long-range pathways and limiting steps.

3. **Insights into Proton Motion:**
   - Flow-Through Centrality results support earlier studies suggesting dual proton motion, with implications for proton−proton correlation.
