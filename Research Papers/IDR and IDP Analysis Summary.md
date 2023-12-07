## Introduction

**Intrinsically Disordered Proteins (IDPs) and Regions (IDRs):** IDPs and IDRs, prevalent in cellular processes, lack a fixed structure and impact signaling and regulation. Very limited resources are available on this domain and exploration can be done.

**Earlier Findings:** Previous work indicated restricted water dynamics around IDPs, emphasizing the role of charged groups. This study delves into the dynamics of hydrogen bonds, considering models for rearrangement and factors like charged group density and solvent-accessible surface area.

**Charged Groups and Surface Area:** 
Concentration of charged groups and the RSA surface area significantly influence hydrogen bond dynamics. 

## Simulation Details

- **Simulation Parameters:**
  - Systems: SUMO-1 and UBC9 simulated with CHARMM 36m force field.
  - Software: GROMACS 2019.6 used for simulations.
  - Initial Configurations: Retrieved from Protein Data Bank (PDB) entries 1A5R and 1A3S.

- **System Preparation:**
  - Solvent: TIP3P water model optimized for CHARMM 36m, with NaCl for 0.15 M salt concentration.
  - Box: Truncated octahedral with 1.4 nm padding from proteins.

- **Analysis of Water Dynamics:**
  - Hydrogen Bond Time Correlation Function: $C(t)$ formula for hydrogen bond existence probability.
  - Lifetime Calculation: τ obtained at$e^{-1}$of$C(0)$via nonlinear least-squares fitting.
  - Average Solvent-Accessible Surface Area (SASA): Computed per residue to assess water access.
  - Retardation Factor: Ratio of$\tau$for polar, cationic, and anionic residues to nonpolar residues.
  - Correlation with Distance: Investigated correlation between hydrogen bond lifetime and residue distance from protein center of mass for small RSA$(<0.3)$.

## Formulas and Definitions

1. **Hydrogen Bond Time Correlation Function (C(t)):**
   -$C(t) = \frac{1}{2} [1 + \text{cos}(\omega(t - t_0))]$, where$\omega = \frac{2\pi}{\tau}$.

2. **Hydrogen Bond Lifetime (τ):**
   -$\tau$is the time where$C(t)$is$e^{-1}$of$C(0)$.

3. **Average Solvent-Accessible Surface Area (SASA):**
   -$\text{RSA}_i = \frac{\text{SASA}_i}{\text{SASA}_{i,\text{max}}}$, where$\text{SASA}_i$is the solvent-accessible surface area for residue$i$and$\text{SASA}_{i,\text{max}}$is the maximum accessible surface area for residue$i$.

4. **Retardation Factor:**
   -$\text{Retardation Factor} = \frac{\tau_{\text{polar}}}{\tau_{\text{nonpolar}}}$,$\frac{\tau_{\text{cationic}}}{\tau_{\text{nonpolar}}}$, or$\frac{\tau_{\text{anionic}}}{\tau_{\text{nonpolar}}}$for polar, cationic, and anionic residues, respectively.

5. **Correlation with Distance:**
   - Correlation between hydrogen bond lifetime and residue distance from the center of mass of the protein for small RSA ($< 0.3$).

These formulas are used to analyze the hydrogen bond dynamics, solvent accessibility, and correlation with protein structure in the context of water interactions with SUMO-1 and UBC9.

## Conclusions

- Study examines hydrogen bond dynamics in proteins (SUMO-1 and UBC9) with intrinsically disordered regions (IDRs) and structured regions.
- Water dynamics around IDPs are constrained, attributed to a higher proportion of charged residues.
- Key factors influencing hydrogen bond dynamics identified:
  - Nature of the residue (anionic, cationic, polar, or nonpolar) impacting water bonding.
  - Time-averaged relative solvent-accessible surface area (RSA) affecting water molecule access.
- Collective effects observed, where a high density of charged residues hinders water motion, leading to longer-lived hydrogen bonds with polar or nonpolar residues.
- Results offer insights into hydrogen bond dynamics on protein surfaces in dilute systems.
- Implications for liquid−liquid phase separation in regions with high densities of IDPs are suggested.
- Further exploration recommended, including study of other SUMO family members and different phases, for a comprehensive understanding.
