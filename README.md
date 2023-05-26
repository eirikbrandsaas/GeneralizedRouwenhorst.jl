# Rouwenhorst

[![Build Status](https://github.com/eirikbrandsaas/Rouwenhorst.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/eirikbrandsaas/Rouwenhorst.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/eirikbrandsaas/Rouwenhorst.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/eirikbrandsaas/Rouwenhorst.jl)

This package implements the generalized Rouwenhorst algorithm from [Fella, G., Gallipoli, G., & Pan, J. (2019). Markov-chain approximations for life-cycle models. Review of Economic Dynamics, 34, 183-201.](https://www.sciencedirect.com/science/article/pii/S1094202519301565?casa_token=1S9HyeijNWoAAAAA:nLK8Xa9hX6bb-uRckgNw8wM6qnVdiEmILVXXuLQlfwp1Ut33Q_-wXm2Bog4MouAzVUqZi3ftHi1x). It is a Julia translation of the [Matlab code](https://github.com/gfell/nsmarkov-matlab) provided by the authors, under the MIT license. If you use this package, I suggest you cite the original authors and their paper.

## To use
```julia
N = 2                       # Nodes
T = 2                       # Time periods
rho = fill(0.95,T)          # Must have length = T
sigma_eps = fill(0.1,T)     # Must have length = T

ygrd, trans = genrouwenhorst(rho,sigma_eps,N,T)
# ygrd is NxT array of grid values, with columns representing time period
# trans is the NxNxT transition matrix. I.e., `trans[1,2,2]' is the probability of transition from node 1 to node 2 in the second period.
```

## Julia code for standard Rouwenhorst methods:
- [QuantEcon.jl](http://quantecon.github.io/QuantEcon.jl/latest/api/QuantEcon.html#QuantEcon.rouwenhorst)
- [Gustavo Pereira](https://github.com/pereiragc/rouwenhorst)