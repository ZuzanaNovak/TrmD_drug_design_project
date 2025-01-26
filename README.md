# TrmD Drug Design Project

This project aims to identify potential inhibitors of TrmD (tRNA (m1G37) methyltransferase), a bacterial enzyme crucial for protein synthesis and a promising antimicrobial drug target. We employ a combined approach of docking simulation and reinforcement learning for ligand generation and optimization.

## Methodology


This project follows a stepwise approach to identify potential TrmD inhibitors by combining reinforcement learning-based molecule generation, docking simulations, and ADMET analysis.

---

### **1. Initial Docking**
- **Objective:** Screen known TrmD-associated ligands to identify suitable candidates for further optimization.  
- **Methodology:**  
  - Docking was performed using AutoDock Vina on 69 ligands from the ChEMBL database.
  - The AdoMet binding pocket of TrmD (structure 4YVG) was used, prepared using AutoDock Tools.
  - Docking results were filtered to select ligands with binding scores ≤ -8 kcal/mol.

---

### **2. Ligand Generation using REINVENT 4**
- **Objective:** Generate and optimize novel molecules based on scaffolds from the initial docking results.  
- **Pipeline:**  
  - **Input:** Selected molecules with good docking scores from the initial docking step.
  - **Stage 1 (Exploration):**
    - Focused on generating diverse molecules using the `mol2mol_scaffold_generic` prior.
    - Employed scoring strategies to optimize molecular weight (200–500 Daltons) and penalize undesirable features.
    - Diversity filters ensured structural novelty.
  - **Stage 2 (Refinement):**
    - Refined molecules from Stage 1, enforcing stricter drug-likeness criteria.
    - Smaller batch sizes and updated diversity filters were used for focused optimization.
  - **Output:** Final library of optimized molecules ready for further docking.

---

### **3. Docking of Generated Molecules**
- **Objective:** Evaluate the binding affinities of REINVENT-generated molecules to TrmD.  
- **Methodology:**  
  - **Initial Screening:**  
    - Docking was performed using AutoDock Vina with moderate exhaustiveness (8).
    - Top 500 molecules (binding scores ≤ -9.052 kcal/mol) were shortlisted.
  - **Refined Docking:**  
    - Shortlisted molecules were re-docked with higher exhaustiveness (32) to refine binding poses.
    - Binding poses were visually inspected to confirm key interactions with active site residues.

---

### **4. ADMET Analysis**
- **Objective:** Assess the pharmacokinetic and safety profiles of the top-ranked ligands.  
- **Tools Used:**  
  - RDKit for calculating basic drug-like properties (e.g., molecular weight, LogP, TPSA).  
  - ADMET-AI, admetSAR, and SwissADME for toxicity predictions (e.g., AMES test, hERG inhibition, hepatotoxicity).  
- **Outcome:**  
  - Molecule 3664 emerged as the top candidate with a strong balance of binding affinity and drug-like properties.
  - Despite favorable ADMET predictions, hepatotoxicity risks were identified, requiring further investigation.

## Results

* **Ligand Generation:** The REINVENT 4 pipeline generated a diverse library of molecules with optimized properties.
* **Docking:**  Several ligands exhibited strong binding affinities to TrmD, indicating potential inhibitory activity.
* **Lead Candidate:** Ligand 3664 emerged as a promising candidate due to its balanced profile of binding affinity, ADMET properties, and structural similarity to a known TrmD inhibitor. However, predicted hepatotoxicity would require further investigation.


## Dependencies

* REINVENT 4
* AutoDock Vina
* RDKit
* ADMET-AI
* admetSAR
* SwissADME
* Open Babel

## Limitations

* ADMET predictions are computational and require experimental validation.
*  Toxicity risks identified for the lead candidate need further investigation and optimization.


## Collaborators
[Ladislav Buček](https://github.com/bucekla)


[Zuzana Nováková](https://github.com/ZuzanaNovak)


[Věra Tereza Štěpánková](https://github.com/stepankverat)
