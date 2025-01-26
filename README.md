# TrmD Drug Design Project

This project aims to identify potential inhibitors of TrmD (tRNA (m1G37) methyltransferase), a bacterial enzyme crucial for protein synthesis and a promising antimicrobial drug target. We employ a combined approach of docking simulation and reinforcement learning for ligand generation and optimization.

## Methodology

The project is divided into two main phases:

**Phase 1: Ligand Generation using REINVENT 4**

* **Data Collection:** Known TrmD-associated ligands are retrieved from the ChEMBL database.
* **REINVENT 4 Pipeline:**  A two-stage process within REINVENT 4 (Mol2Mol mode) is utilized:
    * **Stage 1 (Exploration):** Focuses on generating diverse molecules with general drug-like properties using a scaffold-based prior and multinomial sampling.  Custom alerts and molecular weight optimization guide the generation process.
    * **Stage 2 (Refinement):** Refines the generated molecules by prioritizing those meeting stricter drug-like criteria while maintaining structural diversity. Uses the output of Stage 1 as input.

**Phase 2: Docking Simulation and ADMET Analysis**

* **Initial Docking (Screening):** AutoDock Vina is used for initial screening of the generated ligands against the TrmD binding pocket (identified using P2Rank on the 4YVG structure from *Haemophilus influenzae*). A moderate exhaustiveness parameter is used for faster processing.
* **Final Docking (Refinement):** Top-scoring molecules from the initial docking are re-docked with a higher exhaustiveness parameter for more accurate binding pose prediction.
* **ADMET Analysis:**  The top-ranked molecules undergo ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) analysis using RDKit, ADMET-AI, admetSAR, and SwissADME to assess their drug-like properties and potential risks.

## Results

* **Ligand Generation:** The REINVENT 4 pipeline generated a diverse library of molecules with optimized properties.
* **Docking:**  Several ligands exhibited strong binding affinities to TrmD, indicating potential inhibitory activity.
* **Lead Candidate:** Ligand 3664 emerged as a promising candidate due to its balanced profile of binding affinity, ADMET properties, and structural similarity to a known TrmD inhibitor. However, predicted hepatotoxicity requires further investigation.


## Key Files

* **report.pdf:**  Detailed project report including methodology, results, and discussion.
*  **(Add other relevant files here, e.g., input files for REINVENT 4, docking configuration files, output files containing generated ligands, etc.)**

## Dependencies

* REINVENT 4
* AutoDock Vina
* RDKit
* ADMET-AI
* admetSAR
* SwissADME
* (Add any other dependencies here)


## Usage

**(Provide instructions on how to reproduce the results, if applicable)**


## Limitations

* ADMET predictions are computational and require experimental validation.
*  Toxicity risks identified for the lead candidate need further investigation and optimization.


## Collaborators
