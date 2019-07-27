## [Mesenchymal precursor cells in adult nerves contribute to mammalian tissue repair and regeneration](https://doi.org/10.1016/j.stem.2018.10.024)

Matthew J. Carr, Jeremy S. Toma, Adam P.W. Johnston, Patrick E Steadman, Scott A. Yuzwa, Neemat Mahmud, Paul W. Frankland, David R. Kaplan, Freda D. Miller. Cell Stem Cell, 2019.

Data portal by scClustViz https://baderlab.github.io/scClustViz .

## Abstract

Peripheral innervation plays an important role in regulating tissue repair and regeneration. Here, we provide evidence that injured peripheral nerves provide a reservoir of mesenchymal precursor cells that can directly contribute to murine digit tip regeneration and skin repair. In particular, using single-cell RNA sequencing and lineage tracing we identify transcriptionally-distinct mesenchymal cell populations within the control and injured adult nerve, including neural crest-derived cells in the endoneurium with characteristics of mesenchymal precursor cells. Culture and transplantation studies show that these nerve-derived mesenchymal cells have the potential to differentiate into non-nerve lineages. Moreover, following digit tip amputation, the neural crest-derived nerve mesenchymal cells contribute to the regenerative blastema and ultimately to the regenerated bone. Similarly, neural crest- derived nerve mesenchymal cells contribute to the dermis during skin wound healing. These findings support a model where peripheral nerves directly contribute mesenchymal precursor cells to promote repair and regeneration of injured mammalian tissues.

## Usage

This is the data package to view the single-cell rna-seq data in this paper. Installation (run once) is as follows:
```{r}
install.packages("devtools")
devtools::install_github("millerkaplanlab/MouseSciaticNerve")
```

After installation, the data can be viewed with the viewMouseSciaticNerve() function. For example:
```{r}
library(MouseSciaticNerve)
viewMouseSciaticNerve("Inj9dBeads")
```
