# Fast modular square root

This library implements multi-threaded modular square-root for the HOLMES paper, which is used to prove a Legendre symbol in zero knowledge.

## Requirements

**OpenMP.** Most general-use C/C++ compilers support OpenMP. To install OpenMP for Apple Clang++, use: 

```bash
brew install libomp
```

**FLINT.** Instructions to install the FLINT library can be found [here](https://flintlib.org/doc/building.html).

## Supported prime number field

This library supports the fields of two prime numbers.

The first prime number is `2^62 - 2^16 + 1`, which is the prime used in [our fork of EMP-ZK](https://github.com/holmes-anonymous-submission/emp-zk) for the QuickSilver scheme.

The second prime number is `2^62 - 2^26 - 2^25 + 1`, which is the prime number we prepare for the Virgo scheme.

## Algorithm

We arrange all of the numbers that we wish to find the square roots for as a vector of FMPZ integers. Later, we perform parallel modular square roots over these numbers.

Each thread will perform the following operations to find the square root (mod p). If the number `x` is not a quadratic residue (and does not have a square root), we compute the square root of `QNR * x` where `QNR` is a quadratic nonresidue.
- Assign our number as the variable `x`
- Check if `x` is a quadratic residue using the FLINT Jacobi symbol algorithm
- If the `x` is a quadratic residue, set `A := x`. 
- If the `x` is not a quadratic residue, set `A := x * first found quadratic nonresidue`. 
- Find the square root of `A` using `fmpz_sqrtmod()`.

We find the first quadratic nonresidue by a brute-force search.


## Regulatory issue

This repository is not subject to the U.S. Export Administration Regulation (EAR) because it is publicly available.



For more information about this regulatory issue, see [this post](https://www.eff.org/deeplinks/2019/08/us-export-controls-and-published-encryption-source-code-explained) by Electronic Frontier Foundation (EFF).
