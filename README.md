# Genetic-data-exploration

*Exploring and classifying genetic data.*

In recent years the size and quality of biological data sets have increased at a tremendous rate. It is fair to say that we 
have not yet reached the plateau of that development, so that students and researchers alike will be facing even larger data 
sets in the near future.

As the amount of information increases, so must the methods we use to explore, analyse and even present that data. 
The field of population genetics is undergoing a profound paradigm shift, as researchers move from exploiting the sparse, 
hard-won data sets of old to finding ways to understand and encompass the swaths of data their laboratories turn in. 

In recent years i have had the privilege of working on one of the largest genetic data sets to date (2018). In summary, this 
project involved the description of genetic structure across thousands of possible sites. I would like to leave 
here some of the more useful insight and tools that resulted. This repository consists in a series of Jupyter notebooks, 
each exploring some method or applications i found particularly useful.

Notebooks are organised in order of increasing complexity and, in most cases, of chronological development. This development was driven in 
greater part by the need to study the relation between statistcis that describe structure in data. Metrics to describe structure vary
across fields of research. Here, the focus is on the description of genetic data. 

Use NBviewer directly to explore the [notebook library](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/tree/master/)

## I. Binary samples. Population distribution and identification.

- Notebooks: 1 - 5

Structure in data is either used to study real-world processes influencing the distribution of observed variables or as a basis for
prediction. From a mathematical point of view these two goals are indistinguishable. In either case the focus of scientific research is
on the description of the laws that govern the generation of new data. This description comes in the form of the combined distributions
of all variables available to describe the subject of study. 

In the study of the variation of binary data, the simplest assumption is that of a binomial probability of observation at the smallest 
unit of measurment available. In the case of genetic data, this unit is the single nucleotide polymorphism (SNP).

Througout the notebooks provided, individuals are simulated as samples of *L* binary markers (0, 1, .., i), coded 0 and 1, of frequency *pi* 
within a given population *K*. Populations were simulated as multivariate Bernoulli variables *Sk*, where each variable *Ski* represents the 
binomial probability of an event. We further assumed independent markers and modelled the distributions of allele frequencies *Sk* from the 
Beta distribution, as an approximation of the stationary allele frequencies described in Williamson *et al.* (2004). Allele frequency vectors 
of size *L* were sampled from the Beta distribution for various combinations of mean and variance of this distribution (see [Notebook 1](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/1.%20Generating_haplotypes.ipynb)).

To the extent that individuals are defined as combinations observations of the same weight, principal component analysis (PCA) presents an intuitive
summarisation of population samples. In this context, given a sufficient number of samples, a basic description of a population entity is 
the probability density function of its samples in PCA feature space. In Notebooks [2.](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/2.%20Local_classification.ipynb),
 [3.](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/3.%20Mislabelling.ipynb), [4.](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/4.%20X-material.ipynb),
  and [5.](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/5.%20Visualizing%20KDE.ipynb), the use of kernel
density estimates ([KDE](https://scikit-learn.org/stable/modules/density.html)) for 
the characterisation and assignment of individual samples is explored. 


## II. Population structure.

- Notebooks: 6 - 


[NBviewer](https://nbviewer.jupyter.org/) supports the 2D images produced in the scripts (copy-paste the url of a notebook onto the tab provided).

Not interactive plots however.


