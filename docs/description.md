# Method description

This section describes how the *Risk of non-adaptedness* (RONA) method works, and how it was implemented in *pyRona*.


## Requirements

The RONA method requires a matrix of set grouped individuals characterized for three different components:
* Genetic information (in this case, allele frequencies across a set of loci)
* Current environmental data
* Future predictions of the provided environmental variables


## Overview

After establishing marker-environmental associations, based on **present allele frequencies and present environmental variables**, the inferred linear model is used to estimate **future expected** allele frequencies.
The **average difference** between **current** and **future expected** allele frequencies can be viewed as the required average change in allele frequency for the species to adapt to future local conditions. This value is called the "RONA value".


## Original description

The original RONA method is described in [Rellstab et al. 2016](doi.wiley.com/10.1111/mec.13889). Here is a brief recap of the main points from the original research paper. For a more detailed understanding it is recommended that you actually read the research paper.

In [Rellstab et al. 2016](doi.wiley.com/10.1111/mec.13889), the software *LFMM* [Frichot et al. 2013](doi:10.1093/molbev/mst063)  is used to calculate environmental associations between environmental factors and allele frequencies. In order to reduce the amount of false positives, any associations with a *Spearman's rank coefficient* below 0.3 was removed.

Present and future environmental factors values were then compared using paired *t*-tests and ranked according to the resulting *p*-values. Only the **three** most differentiated factors were further analysed.

To estimate the adaptedness of populations to future local conditions, the extent to which the present allele frequencies at climate-related SNPs differ in average from those expected under modelled future climate, given the respective linear model of environmental association was assessed. This required average change in allele frequency is defined as the risk of nonadaptedness (RONA) to future conditions.

To infer RONA, first, simple linear regressions of the selected environmental factors of the present climate (EF present) with the alternative allele frequencies (compared to the reference sequence, AAF present) for all loci were performed.

<!-- Per environmental factor tested, we then selected, if possible, the top 20 loci from significant (P < 0.05) linear regressions that were also candidate loci in the LFMM analysis. These loci are the same for all three species. To exclude linked loci, only one locus per target was included. -->
