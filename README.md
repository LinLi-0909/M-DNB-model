#  M-DNB Factors Orchestrate Cell Fate Determination at Tipping Points during Mesendodermal Differentiation of Human Embryonic Stem Cells 
## Overview
## Contents

- [Overview](#overview)
- [M-DNB analysis Guide](./LICENSE)
- [License](./LICENSE)
- [Citation](#citation)
- [Contact](#Contact)

## Overview
The generation of ectoderm, mesoderm, and endoderm layers is the most critical biological process during the gastrulation of embryo development. Such a differentiation process in human embryonic stem cells (hESCs) is an inherently nonlinear multi-stage dynamical process which contain multiple tipping points playing crucial roles in the cell-fate decision. However, the tipping points of the process are largely unknown, letting alone the understanding of the molecular regulation on these critical events. Here by designing a module-based dynamic network biomarker (M-DNB) model, we quantitatively pinpointed two tipping points of the differentiation of hESCs towards definitive endoderm, which leads to the identification of M-DNB factors (FOS, HSF1, MYCN, TP53 and MYC) of this process. 
The stage-specific and essential roles of M-DNB factors in the cell-fate decision were confirmed by the differentiation experiments. We demonstrate that before the tipping points, M-DNB factors are able to maintain the cell states and orchestrate cell fate determination during hESC (ES)-to-ME and ME-to-DE differentiation processes, which not only leads to better understanding of endodermal specification of hESCs but also reveals the power of the M-DNB model to identify critical transition points with their key factors in diverse biological processes including cell differentiation and transdifferentiation dynamics.

## M-DNB analysis Guide
### 1. Get network from PPI network
First, we should obtain the network from PPI network from STRING (https://cn.string-db.org/).
The PPI network of human can be downloaded from https://github.com/LinLi-0909/M-DNB-model/blob/main/9606.protein.links.v10.zip.
### 2. Identify critical points and CI values of each gene module
**CI = Get_CI(data,time_Idx,feature,network);** <br />
**data is time-series scRNA-seq dataset;** <br /> genes*cells, e.g. Chu-time dataset(download from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE75748&format=file&file=GSE75748%5Fsc%5Ftime%5Fcourse%5Fec%2Ecsv%2Egz).<br />
![image](https://user-images.githubusercontent.com/63344240/200508091-bc34407b-5bfa-4942-bfa6-e53e28df450d.png)
**time_Idx is time points of samples;** <br /> e.g Chu-time dataset contains 6 time points.The mat file can be downloaded from https://github.com/LinLi-0909/M-DNB-model/blob/main/timeIdx.mat
![image](https://user-images.githubusercontent.com/63344240/200511593-cda817ff-8ada-4033-9839-0b803e68e7f3.png)

**feature is genes in dataset.** <br /> e.g.row name of Chu-time dataset.
**network is PPI network obtained from STRING.** 

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
