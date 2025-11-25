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

## System Requirements & Dependencies

The analysis was performed using **MATLAB** and **Python**. To reproduce the results, the following dependencies are required:

### MATLAB
* **MATLAB Version:** (e.g., R2021b or later)
* **Toolboxes:** * SPM12
    * DPABI / DPARSF
    * Brain Connectivity Toolbox (BCT)
    * (Please list any other specific MATLAB toolboxes used here)

### Python
* **Python Version:** (e.g., 3.8+)
* **Libraries:**
    * `numpy`
    * `pandas`
    * `nibabel`
    * `scipy`
    * (Please list any other specific Python libraries used here)

---

## Usage

This repository is organized to facilitate the reproduction of the analysis pipeline described in the manuscript.

1.  **Preprocessing:** Scripts located in the `preprocessing` folder handle the initial fMRI data processing.
2.  **Energy Analysis:** The core analysis based on network control theory is located in the `analysis` folder.
3.  **Visualization:** Scripts used to generate the figures in the main text are provided in the `plotting` folder.

---

## Citation

If you use this code in your research, please cite our paper:

> Liu, Q., Xiong, H., et al. (2025). Energy inefficiency underpinning brain state dysregulation in individuals with major depressive disorder. *Nature Mental Health*. (In Press)
