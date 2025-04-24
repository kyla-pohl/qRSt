# qRSt
This is my implementation of a generalized version of the Robinson-Schensted algorithm as described by Aigner and Frieden in their paper [qRSt: A probabilistic Robinson–Schensted correspondence for Macdonald polynomials](https://arxiv.org/pdf/2104.13846). The implementation is provided in a Jupyter notebook with a Sage kernel.

The celebrated Robinson-Schensted algorithm provides an explicit bijection from permutations to pairs of standard Young tableaux of the same shape. In their 2021 paper, Aigner and Friedan describe a $q,t$-generalization of RS which they call $qRSt$. It specializes to several known RS generalizations, in particular column-insertion RS (when $q \to \infty$, $t \to 0$) and row-insertion RS (when $q, t \to 0$). 

In this repository, I have implemented $qRSt$ in SageMath. The qRSt_functions.sage file contains all of the functions needed for the algorithm and is safely ignored. Put in your desired symmetric group element and parameters $q$ and $t$ in the qRSt.ipynb file and run the notebook. The output will be two probabilistically determined standard young tableaux of the same shape.
