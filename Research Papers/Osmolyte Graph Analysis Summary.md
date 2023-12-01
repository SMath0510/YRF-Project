## Simulation Details

- **Software:** GROMACS 2016.4
- **Systems Studied:** Neat water and osmolyte solutions (ethanol, glycerol, α-glucose, β-glucose, sorbitol, trehalose, urea)
- **Osmolyte Concentration:** 20 wt%
- **Force Field:** CHARMM36
- **Water Models:** SPC/E and CHARMM-modified TIP3P

## Ratios Defined for Network Analysis

1. **Coordination Number**
   - Number of nearest neighbors.

2. **Preferential Contact Ratio**
   - Local Osmolyte to water coordination number normalized by the total number of osmolyte atoms that can form hydrogen bonds to the total number of water molecules.

## Simulation Parameters

- **Cutoff for Electrostatic Interactions:** 1.2nm.
- **Hydrogen Bond Conditions:**
  - Donor Acceptor angle < 30.
  - Weight of the bond depends on the angle and the distance (assumed to be 1 in the graph analysis point of view).

## Network Analysis Metrics

We calculate the following centrality metrics for the graphs:

1. **Betweenness Centrality:**
   - A measure of the number of shortest paths that pass through a node in the graph. Nodes with high betweenness centrality play a crucial role in maintaining connectivity within the network.

2. **Closeness Centrality:**
   - A measure of how quickly a node can access other nodes in the network. Nodes with high closeness centrality are close to many other nodes in the graph.


## Formulas and Definitions

### Coordination Number and Preferential Contact Ratio Calculation:

#### Coordination Number ($n_i$):
$ n_i = 4 \pi \rho_i,b \int_{r=0}^{R} g_{wi}(r) r^2 \, dr $

- $i$ corresponds to either water ($wat$) or osmolyte ($osm$),
- $R$ is the cutoff distance of 0.35 nm.

#### Bulk Density (ρ):
$ \rho = \frac{N}{V} $

- $N$ is the total number of molecules,
- $V$ is the simulation box volume.

#### Preferential Contact Ratio ($γ$):
$ γ = \frac{n_{osm}}{n_{wat}} \frac{N_{wat} - 1}{N_{osm} - Nh} $

- $N_h$ is the number of heavy atoms per osmolyte.

### Interaction Energy Calculation:

$ Interaction \ Energy = Dispersion \ Interaction + Electrostatic \ Interaction $

### Tetrahedral Order Parameter (q4,i):

$ q_{4,i} = \frac{1}{8} \sum_{j,k} \cos^2 \theta_{jik} $

### Hydrogen Bond Strength Calculation:

$ w = \cos^2 \theta \left(1 - \frac{r - r_{min}}{r_{max} - r_{min}}\right) $

- $r$ is the distance between donor and acceptor atoms,
- $θ$ is the angle between the hydrogen atom, the donor, and the acceptor.

### Graph Network Analyses:

#### Closeness Centrality ($Cc$):

$ Cc(u) = \frac{n}{N(N-1)} \sum_{v} \frac{1}{d(u,v)} $

- $n$ is the number of reachable nodes to the node $u$,
- $d(u,v)$ is the shortest distance between the node $u$ and a reachable node $v$.

#### Betweenness Centrality ($Cb$):

$ Cb(u) = \frac{\sum_{s,t} np(s,t|u)}{(N-1)(N-2)} $

- $n_p(s,t)$ is the number of shortest paths between nodes $s$ and $t$,
- $n_p(s,t|u)$ is the number of such paths that pass through the node $u$.

### Network Properties Comparison:

#### Relative Mean Path Length ($Δp$):
$ Δp = \langle p_{soln} \rangle - \langle p_{NW} \rangle $

#### Effective Betweenness/Closeness Centrality ($η$):
$ η = \frac{\langle Cb_{osm} \rangle}{\langle Cb_{wat} \rangle} $

#### Effective Closeness ($ξ$):
$ ξ = \frac{\langle Cc_{wat} \rangle - \langle C̃c_{wat} \rangle}{\langle Cc_{wat} \rangle} $

## Conclusions

1. **Water Structure Perturbation:** Binary solutions show slight changes in water structure, affecting water−water hydrogen bonds and tetrahedral order parameters, but these are compensated by water−osmolyte hydrogen bonds.

2. **Osmolyte Interaction Preference:** Osmolytes preferentially hydrogen bond with water, replacing some water molecules while retaining the local structure.

3. **Thermodynamic Stability:** Thermodynamic properties of water in binary solutions remain unaltered, with osmolyte-water interactions crucial for solution properties.

4. **Graph Theory Analysis:** Graph theory reveals binary solution hydrogen-bond networks are percolated, exhibiting shorter mean paths compared to neat water.

5. **Osmolyte Connectivity:** Osmolytes act as hubs in the solution hydrogen-bond network, influencing graph properties based on their hydrogen-bonding capacity.

6. **Effect of Temperature:** Minimal temperature impact on graph properties, but reduced translational motion stiffens the hydrogen-bond network.

7. **Urea and Glycerol Similarities:** Similar local and graph properties in urea and glycerol solutions suggest classifying osmolytes as denaturant or stabilizer is challenging based solely on hydrogen-bond analysis.

8. **Direct Protein-Osmolyte Interaction:** The observed correlation between osmolyte binding energy, their hydrogen bonds, and protein conformational changes suggests direct protein-osmolyte interactions drive structural alterations.

9. **Hydrophilic Interactions Significance:** The study emphasizes the role of hydrophilic interactions in protein folding, with network-based analysis providing insights into osmolyte influence on protein structure.

10. **Protein Conformational Equilibria:** The extent of hydrogen bonding between proteins and osmolytes could influence protein conformational equilibria, contributing to the understanding of protein folding dynamics.