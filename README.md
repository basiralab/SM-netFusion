# SM-netFusion (Supervised multi-topology network cross-diffusion) 

SM-netFusion (Supervised multi-topology network cross-diffusion) for a fast and accurate graph data classification code, created by Islem Mhiri. Please contact islemmhiri1993@gmail.com for inquiries. Thanks.

Although limited, existing brain network atlases (BNA) rely on a similarity network diffusion and fusion technique, which only considers node degree as a topological measure in the cross-network diffusion process, thereby overlooking rich topological measures of the brain network (e.g., centrality). Moreover, both diffusion and fusion techniques are implemented in fully unsupervised manner, which might decrease the discriminative power of the estimated BNAs. To address these issues, we design a simple but effective  a supervised multi-topology network cross-diffusion (SM-netFusion) framework for learning a BNA which satisﬁes the following constraints:(i) it is well-representative that consistently captures the unique and distinctive traits of a population of functional networks, (ii) it is well-centered that occupies a center position optimally near to all individuals, and (iii) it reliably identiﬁes the most discriminative disordered brain connections by comparing templates estimated using disordered and healthy brains, respectively.


![SM-netFusion pipeline](https://github.com/basiralab/SM-netFusion/blob/master/Fig1.png)

# Detailed proposed SM-netFusion pipeline

This work will be published in **MICCAI 2020 LNCS Springer Proceedings** (https://www.miccai2020.org/en/).

**Supervised multi-topology network cross-diffusion (SM-netFusion)** presents the ﬁrst work for supervised network cross-diffusion based on graph topological measures (SM-netFusion) by enhancing the non-linear fusion process using a weighted mixture of multi-topological measures.  Our learning-based framework comprises three key steps. (1) Class-speciﬁc feature extraction and clustering, (2) Class-speciﬁc supervised multi-topology network cross-diffusion, (3) Identiﬁcation of the discriminative connectional ﬁngerprint. Experimental results and comparisons with the state-of-the-art methods demonstrate that SM-netFusion can achieve the best results in terms of centerdness, representativness and discriminativness and further boosted classification accuracy. We evaluated our proposed framework from ABIDE preprocessed dataset (http://preprocessed-connectomes-project.org/abide/).


# Demo
The code has been tested with MATLAB 2020a on Windows 10. GPU is not needed to run the code.
In this repository, we release the SM-netFusion source code trained and tested on a simulated heterogeneous graph data from 2 Gaussian distributions as shown below:

**Data preparation**

We simulated random graph dataset from two Gaussian distributions using the function simulateData.m. The number of graphs in class 1, the number graphs in class 2, and the number of nodes (must be >20) are manually inputted by the user when starting the demo.

To train and evaluate SM-netFusion code on other datasets, you need to provide:
• A tensor of size (((n/5) × 4 )× m × m) stacking the symmetric matrices of the training subjects. n denotes the total number of subjects and m denotes the number of nodes.<br/>
• A vector of size ((n/5) × 4 )stacking the training labels.<br/>
• The number of selected features Nf.<br/>
• A boolean variables ‘displayResults’ ∈ [0, 1].<br/>

If displayResults = 1 ==> display (Atlas of group 1, Atlas of group 2, top features matrix and the circular graph at each cross-validation run). This is cool but it might slow down the demo a bit. <br/>
If displayResults = 0 ==> no display except for the average results across all cross-validation runs.

**The SM-netFusion outputs are:**
• A matrix of size (m × m) storing the network atlas of group 1.<br/>
• A matrix of size (m × m) storing the network atlas of group 2.<br/>
• A vector of size (Nf × 1) stacking the indices of the top discriminative features.<br/>

**Train and test SM-netFusion**

To evaluate our framework, we used 5-fold cross validation strategy.

To try our code, you can use: SM-netFusion_demo.m

# Example results

Set the boolean variable displayResults=1 and check our results. You can set the number of samples (i.e., graphs) from class 1 to 30, from class 2 to 30, and the size of each graph to 60 (nodes).


# YouTube videos on SM-netFusion

To install and run SM-netFusion, check the following YouTube video:

https://www.youtube.com/watch?v=-dDn8CT4mH0

To learn about how SM-netFusion works, check the following YouTube video:

https://www.youtube.com/watch?v=eWz65SyR-eM

# Acknowledgement

We used the following codes:

EasyMKL code from https://github.com/okbalefthanded/EasyMKL 

SIMLR code from https://github.com/BatzoglouLabSU/SIMLR/tree/SIMLR/MATLAB

SNF code from http://compbio.cs.toronto.edu/SNF/SNF/Software.html

CircularGraph code from https://www.github.com/paul-kassebaummathworks/circularGraph



# Related references

Easymkl: a scalable multiple kernel learning algorithm: Aiolli, F., Donini, M. Neurocomputing169(2015) 215–224. [https://www.sciencedirect.com/science/article/abs/pii/S0925231215003653]

Similarity Network Fusion (SNF): Wang, B., Mezlini, A.M., Demir, F., Fiume, M., Tu, Z., Brudno, M., HaibeKains, B., Goldenberg, A., 2014. Similarity network fusion for aggregating data types on a genomic scale. [http://www.cogsci.ucsd.edu/media/publications/nmeth.2810.pdf] (2014) [https://github.com/maxconway/SNFtool].

Single‐cell Interpretation via Multi‐kernel LeaRning (SIMLR): Wang, B., Ramazzotti, D., De Sano, L., Zhu, J., Pierson, E., Batzoglou, S.: SIMLR: a tool for large-scale single-cell analysis by multi-kernel learning. [https://www.biorxiv.org/content/10.1101/052225v3] (2017) [https://github.com/bowang87/SIMLR_PY].

Paul Kassebaum (2019). circularGraph (https://www.github.com/paul-kassebaum-mathworks/circularGraph), GitHub. Retrieved December 26, 2019

# Please cite the following paper when using SM-netFusion:

@article{mhiriMICCAI2020,
  title={Supervised Multi-topology Network Cross-diffusion for Population-driven Brain Network Atlas Estimation},<br/>
  author={Mhiri, Islem, Ben Khelifa, Anouar, Mahjoub, Mohamed Ali and Rekik, Islem},<br/>
  booktitle={International Conference on Medical Image Computing and Computer-Assisted Intervention},<br/>
  pages={},<br/>
  year={2020},<br/>
  organization={Springer}<br/>
}<br/>

# License
Our code is released under MIT License (see LICENSE file for details).

# Contributing
We always welcome contributions to help improve NAG-FS and evaluate our framework on other types of graph data. If you would like to contribute, please contact islemmhiri1993@gmail.com. Many thanks.
