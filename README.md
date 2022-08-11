# mnempy

## Files:

main.py contains the main nem and mnem functions and a function to simulate data

low.py contains internal helper functions

## Install:

Install within its directory.

```
pip install -e .
```

## Usage:

Create data from two networks with three gene, two effectors per gene and 15 samples.

```
import mnempy as mn

n = 3
K = 2
m = 2
s = 15

sim = mn.sim(K, n, m, s)
```

We extract the data object and the perturbation map (which gene has been perturbed
in which sample) from the simulation object.

```
D = sim["D"]
rho = sim["rho"]
```

Finally, we use the mnem algorithm to compute the mixture and can compare to the
ground truth.

```
result = mn.mnem(D, rho, K)
result["phis"]
sim["phis"]
```

Several runs are suggested, depending on noise and complexity of the ground truth. If
the number of components K is unknown, several different Ks are suggested. However, the
log likelihood cannot be compared directly. Some form of Information Criterion like BIC
or AIC (not implemented) is necessary to select the best K.

DISCLAIMER:
-----------

This program is fully functional and computes a mixture of directed graphs based on
the M&NEM algorithm (Pirkl & Beerenwinkel. 2018).

However, this program is written in python as a programming exercise and should
only be used for testing purposes (feedback is appreciated). For real data analysis
I recommend the more thoroughly tested R/Bioconductor package mnem (https://github.com/cbg-ethz/mnem).

## References:

Pirkl, M., Beerenwinkel, N.; Single cell network analysis with a mixture
of Nested Effects Models, Bioinformatics, Volume 34, Issue 17, 1 September
2018,
Pages i964-i971, https://doi.org/10.1093/bioinformatics/bty602.