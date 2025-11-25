# Energy inefficiency underpinning brain state dysregulation in individuals with major depressive disorder

This repository contains the custom code and analysis scripts for the manuscript **"Energy inefficiency underpinning brain state dysregulation in individuals with major depressive disorder"**.

## Authors
Qianhui Liu, Hui Xiong, Weiyang Shi, Shiqi Di, Xinle Cheng, Nianyi Liu, Na Luo, Yu Zhang, and Tianzi Jiang

**Correspondence:** Prof. Tianzi Jiang (jiangtz@nlpr.ia.ac.cn)  
*Brainnetome Center, Institute of Automation, Chinese Academy of Sciences*

---

## Abstract

Disruptions in brain state dynamics are a hallmark of major depressive disorder (MDD), yet their underlying mechanisms remain unclear. Building on network control theory, this study revealed that energy inefficiency, characterized by elevated energy costs and reduced control stability, drove decreased state stability and increased state-switching frequency in MDD. 

Key brain regions, including left dorsolateral prefrontal cortex and insula, exhibited impaired energy regulation capacity (a metric validated against cerebral metabolism). Moreover, these region-specific energy patterns were correlated with depressive symptom severity. Neurotransmitter and gene expression association analyses linked these energy deficits to intrinsic biological factors, notably the 5-HT2a receptor and astrocytes. 

These findings shed light on the energetic mechanism underlying brain state dysregulation in MDD and its associated biological underpinnings, highlighting brain energy dynamics as a potential biomarker by which to explore therapeutic targets and advance precise interventions for restoring healthy brain dynamics in depression.

---

## Repository Structure

The code is organized into three main modules corresponding to the analytical steps in the study:

### 1. [`Data_Preprocessing/`](Data_Preprocessing/)
* **Purpose:** Preprocessing of fMRI/dMRI data and construction of the Brainnetome Atlas (BNA) based networks.
* **Key Scripts:** Pipeline adapted from UKB-connectomics for unified cortical-subcortical parcellation.

### 2. [`Brain_State_Analysis/`](Brain_State_Analysis/)
* **Purpose:** Extraction of recurrent brain states and calculation of temporal dynamic metrics.
* **Key Metrics:** Fractional Occupancy (FO), Dwell Time (DT), Appearance Rate (AR), Transition Probability (TP).
* **Includes:** Validation scripts for "NoLIM" (excluding Limbic network) conditions.

### 3. [`Energy_Analysis/`](Energy_Analysis/)
* **Purpose:** Quantification of control energy using Network Control Theory (NCT).
* **Key Metrics:** Global Transition Energy, Regional Energy Regulation Capacity (rERC).
---

## External Resources & Dependencies

This study relies on several open-source libraries. We gratefully acknowledge the authors of the following tools:
* **Preprocessing:** [UKB-connectomics](https://github.com/sina-mansour/UKB-connectomics), FSL, Freesurfer, MRtrix3.
* **Brain States:** [control_costs](https://github.com/NeuroenergeticsLab/control_costs), [energy_landscape](https://github.com/singlesp/energy_landscape).
* **Control Energy:** [nctpy](https://github.com/LindenParkesLab/nctpy).
* **Statistical Analysis:** * PLS Correlation: [myPLS](https://github.com/danizoeller/myPLS) (modified for spin tests).
    * PLS Regression: [NSPN_WhitakerVertes_PNAS2016](https://github.com/KirstieJane/NSPN_WhitakerVertes_PNAS2016).
    * Spin Tests: [rotate_parcellation](https://github.com/frantisekvasa/rotate_parcellation), [ENIGMA Toolbox](https://github.com/MICA-MNI/ENIGMA).
    * Meta-analysis: [NiMARE](https://github.com/neurostuff/NiMARE).
* **Gene & Cell Analysis:** [Metascape](https://metascape.org), [hansen_genescognition](https://github.com/netneurolab/hansen_genescognition).
* **Null Models:** [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet), [Geometry-preserving nulls](https://www.brainnetworkslab.com/coderesources).
* **Visualization:** [Neuromaps](https://github.com/netneurolab/neuromaps).

## Citation

If you use this code or the associated findings, please cite:

> Liu, Q., Xiong, H., et al. (2025). Energy inefficiency underpinning brain state dysregulation in individuals with major depressive disorder. *Nature Mental Health*. (In Press)
