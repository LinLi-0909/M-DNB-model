#  M-DNB Factors Orchestrate Cell Fate Determination at Tipping Points during Mesendodermal Differentiation of Human Embryonic Stem Cells 
## Overview
## Contents

- [Overview](#overview)
- [M-DNB analysis Guide](#M-DNB analysis Guide)
- [License](./LICENSE)
- [Citation](#citation)
- [Contact](#Contact)

## Overview
The generation of ectoderm, mesoderm, and endoderm layers is the most critical biological process during the gastrulation of embryo development. Such a differentiation process in human embryonic stem cells (hESCs) is an inherently nonlinear multi-stage dynamical process which contain multiple tipping points playing crucial roles in the cell-fate decision. However, the tipping points of the process are largely unknown, letting alone the understanding of the molecular regulation on these critical events. Here by designing a module-based dynamic network biomarker (M-DNB) model, we quantitatively pinpointed two tipping points of the differentiation of hESCs towards definitive endoderm, which leads to the identification of M-DNB factors (FOS, HSF1, MYCN, TP53 and MYC) of this process. 
The stage-specific and essential roles of M-DNB factors in the cell-fate decision were confirmed by the differentiation experiments. We demonstrate that before the tipping points, M-DNB factors are able to maintain the cell states and orchestrate cell fate determination during hESC (ES)-to-ME and ME-to-DE differentiation processes, which not only leads to better understanding of endodermal specification of hESCs but also reveals the power of the M-DNB model to identify critical transition points with their key factors in diverse biological processes including cell differentiation and transdifferentiation dynamics.

## M-DNB analysis Guide
### 1. Get network from PPI network
First, we should obtain the network from PPI network from STRING (https://cn.string-db.org/)
### 2. Identify critical points and CI values of each gene module

CI = Get_CI(data,time_Idx,feature,whole);
data is time-series dataset
time_Idx is time points of samples
feature is gene considered in dataset
whole is PPI network obtained from STRING

### 3. Identify DNB genes at critical points
 [topmCI,QI]=Get_Critical_Indicators(timeIdx,CI,m);
 m is the top m DNB genes with CI in given critical point.
### 4. Find M-DNB fatcors (TFs of DNB genes)
We can find TFs of DNB genes based on IPA.

## Contact
Please contact us:  
Lin Li: lilin6@sibcb.ac.cn

## Copyright
©2022 Lin Li [Chen Lab]. All rights reserved.
