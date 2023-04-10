<!--
Copyright (C) 2008-2014,2023 Mitsubishi Electric Research Laboratories (MERL)

SPDX-License-Identifier: AGPL-3.0-or-later
-->
# Parallel Quadratic Programming

## Summary

This is a baseline implementation of the Parallel Quadratic Programming (PQP) algorithm, which in its simplest form is a multiplicative fixpoint for the non-negative vector x>=0 that minimizes the quadratic form `f(x)==x'*Q*x/2-h'*x`, where Q is symmetric semi-definite and h is any finite real vector.  The fixpoint works by splitting f(x) into a difference of two quadratic forms, both of which are guaranteed to have strictly non-negative gradients. The elementwise ratio of these gradients is used to scale the solution estimate x, also elementwise; so all elements of x can be updated in parallel, synchronously or asynchronously, via very simple operations: two inner products, one scalar division, and one scalar multiply.   This multiplicative update is guaranteed to improve the value of f(x).  For positive definite Q, repeated updates will converge at a linear rate to the optimal x from any strictly positive initial guess. For semidefinite Q, the rate of convergence will depend on the splitting. Three of many possible splittings are implemented below and can be chosen by uncommenting lines of code.  The simplest, maximally sparse split yields a one-line algorithm:

 iterate: `x = x.*(max(-Q,0)*x+max(h,0))./(max(Q,0)*x+max(-h,0));`

This code can also be used to solve general inequality-constrained quadratic programs in their KKT dual form.  The algorithm can be quite fast and has been used solve very large problems at real-time rates for active computer vision, radiation therapy planning and adaptive delivery, and optimal control of machines including roadway vehicles.

## Installation

As of 2010, this standalone script works in both matlab and octave.  Put it in the current working directory or any directory in the scripts path.

## Usage

Type `help pqp_simple.m` for usage.

## Citation

If you use the software, please cite the following [publications](http://www.merl.com/publications/TR2011-064):

```
@inproceedings{Brand2011sep,
author = {Brand, M. and Chen, D.},
title = {Parallel Quadratic Programming for Image Processing},
booktitle = {IEEE International Conference on Image Processing (ICIP)},
year = 2011,
pages = {2261--2264},
month = sep,
doi = {10.1109/ICIP.2011.6116089},
url = {http://www.merl.com/publications/TR2011-064}
}

@article{DiCairano2013jul,
author = {{Di Cairano}, S. and Brand, M. and Bortoff, S.A.},
title = {Projection-free Parallel Quadratic Programming 
         for Linear Model predictive Control},
journal = {International Journal of Control},
year = 2013,
month = jul,
url = {http://www.merl.com/publications/TR2013-059}
}

@inproceedings{DiCairano2013dec2,
author = {{Di Cairano}, S. and Brand, M.},
title = {On a Multiplicative Update Dual Optimization Algorithm
         for Constrained Linear MPC},
booktitle = {IEEE Conference on Decision and Control (CDC)},
year = 2013,
month = dec,
url = {http://www.merl.com/publications/TR2013-108}
}
```

## Contact

[Matt Brand](http://www.merl.com/people/brand)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for our policy on contributions.

## License

Copyright (c) 2008-2014,2023 Mitsubishi Electric Research Laboratories (MERL).
