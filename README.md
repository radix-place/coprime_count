# Coprime decomposition counts and diagrams

This repository contains the reference implementation accompanying the article:

> *A Functional Framework for Exact Counting of Restricted Coprime Representations of Even Integers*  
> Andrés M. Salazar  
> Pontificia Universidad Javeriana Cali — Colombia  
> (preprint, arXiv)

The code provides:

- An explicit closed formula for the function $g(2n,p)$ (see Theorem 4.13 in the paper),
- A brute-force oracle for direct verification,
- Systematic computational validation,
- Visualization scripts for generating coprime decomposition diagrams.

---

## Contents

- `coprime_count.py`  
  Core implementation of the closed formula for `g(2n,p)`, together with a strict brute-force oracle and pointwise verification routines.

- `prime_test.py`  
  Efficient primality testing based on a precomputed sieve and memory-mapped access.

- `sieve_creator.py`  
  Script to generate the prime sieve file `sieve.npy`.  
  *(Run this script first to create `sieve.npy`, which is required by the other scripts.)*

- `test_count.py`  
  Exhaustive pointwise comparison between the closed formula and the brute-force oracle.

- `graph.py`  
  Script to generate coprime decomposition diagrams.

---

## Requirements

- Python ≥ 3.9  
- NumPy  
- Matplotlib  
