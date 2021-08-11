Usage of the random forrest classifier to predict ovarian cancer subtype
of new sample
================
Bethany M Barnes
01/07/2021

## Background

Epithelial ovarian cancer (OC) is a heterogenous disease consisting of
five major histologically distinct subtypes: high-grade serous (HGSOC),
low-grade serous (LGSOC), endometrioid (ENOC), clear cell (CCOC) and
mucinous (MOC). Although HGSOC is the most prevalent subtype,
representing 70–80% of cases, a 2013 landmark study by Domcke et al.,
found that the most frequently used OC cell lines are not molecularly
representative of this subtype. This raises the question, if not HGSOC,
from which subtype do these cell lines derive? Indeed, non-HGSOC
subtypes often respond poorly to chemotherapy; therefore, representative
models are imperative for developing new targeted therapeutics.

## Methods

Non-negative matrix factorisation (NMF) was applied to transcriptomic
data from 44 OC cell lines in the Cancer Cell Line Encyclopedia,
assessing the quality of clustering into 2–10 groups. Epithelial OC
subtypes were assigned to cell lines optimally clustered into five
transcriptionally distinct classes, confirmed by integration with
subtype-specific mutations. A transcriptional subtype classifier was
then developed by trialling three machine learning algorithms using
subtype-specific metagenes defined by NMF. The ability of classifiers to
predict subtype was tested using RNA-sequencing of a living biobank of
patient-derived OC models. Results Application of NMF optimally
clustered the 44 cell lines into five transcriptionally distinct groups.
Close inspection of orthogonal data sets revealed this five-cluster
delineation corresponds to the five major OC subtypes. This NMF-based
classification validates the Domcke et al. analysis, in identifying
lines most representative of HGSOC, and additionally identifies models
representing the four other subtypes. However, NMF of the cell lines
into two clusters did not align with the dualistic model of OC and
suggests this classification is an oversimplification. Subtype
designation of patient-derived models by a Random Forest transcriptional
classifier aligned with prior diagnosis in 76% of unambiguous cases. In
cases where there was disagreement, this often indicated potential
alternative diagnosis, supported by review of histological, molecular
and clinical features.

## Conclusion

This robust classification informs selection of the most appropriate
models for all five histotypes. Following further refinement on larger
training cohorts, transcriptional classification may represent a useful
tool to support classification of new model systems of OC subtypes.

## Usage

An example R script showing the classification of the cell line
KURAMOCHI is provided. To test a novel cell line, read counts should be
substituted with thhose of the cell line of interest.

The classifier was trained using 44 epithelial ovarian cancer cell lines
with gene annotations from gencode v32. Therefore, the input gene matrix
would ideally have been pre-processed in the same manner, with
normalised reads outputted from e.g. VST function of deseq2.
