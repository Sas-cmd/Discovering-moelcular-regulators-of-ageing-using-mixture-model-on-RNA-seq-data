
**The application of mixture models to RNA-seq data to discover ageing regulators** 
=========================



**Motivation/Abstract**
-----------

Ageing is a complex process. The combined effects of environmental and genetic factors make it challenging to isolate specific regulators. Given the dynamic nature of gene expression, a gene expression can follow different distributions during the aging process. 


![Distribution_SFig1](https://user-images.githubusercontent.com/52276989/156717081-56d82e51-4087-4a0e-a8e6-b68f4b500d19.png)


We can capture the biological variability using mixture models. This is done by modelling the variability via multimodality using multiple different distributions at the gene level for RNA-sequencing (RNA-seq) data. 

We used the [Genotype-Tissue Expression (GTEx)](https://gtexportal.org/home/) cohort to identify lists of candidate genes that clustered according to multimodal distributions with donors that showed significant changes in age. MTOR was the only age-related gene that was identified through our mixture model analysis and not captured through differential expression analysis. We identified mixture model only genes that were common across different tissues, suggesting the presence of systemic ageing genes.  Gene set over representation using the mixture model only genes and the standard differentially expressed gene list resulted in similar pathways, indicating that mixture models are detecting different genes in the same pathway.
The results indicate that modelling gene expression variability using mixture models in conjunction with standard differential gene expression can help uncover regulators that have a potential role in understanding human ageing.


![Figure_1_Flowchart](https://user-images.githubusercontent.com/52276989/156529152-e217b5eb-2c84-4380-8b89-2b08c51eaaf4.svg)

About Data
-----------


Due to the resources required to run mixture models on the GTEx datasets we have provided the clustering outputs from the best mixture models as well as data fom the clustering outputs from all mixture model packages we have run on Github. If data is required and not present please do not hesistate to contact me.

We have also provided the code to download and normalise the GTEx data we used from the [YARN R package](https://github.com/QuackenbushLab/yarn).
