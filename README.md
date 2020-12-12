# GI_InteractionTests
Genetic Interaction tests used to detect mutually exclusive and co-occurring altered gene pairs in childhood cancers.

For more details about both tests we refer to the accompanying manuscript **A comprehensive map of genetic interactions in childhood cancer reveals multiple underlying biological mechanisms** (https://doi.org/10.1101/2020.11.17.385120)

In the **wesme** and **permutation** directory you can find manuals how to run the tests on a local machine or on a hpc cluster.

**WeSME test**

The original WeSME test is published in **WeSME: uncovering mutual exclusivity of cancer drivers and beyond** by Kim et al. in 2017 (https://doi.org/10.1093/bioinformatics/btw242) and the original code can be found on
https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/index.cgi#wesme. 

We have adapted parts of the code for our project in order to run PAN cancer tests. We also changed the way p-values and FDR estimates were calculated. For example, originally only a part of the p-values (those below a certain threshold) were included for FDR estimation. We now include all p-values.

**Permutation test**

The permutation test is based on the genetic interaction test described in **Cancer type-dependent genetic interactions between cancer driver alterations indicate plasticity of epistasis across cell types** by Park & Lehner in 2015 (https://doi.org/10.15252/msb.20156102).
