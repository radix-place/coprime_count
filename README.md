# coprime_count.py

**Exact functional formulas for restricted coprime representations of even integers**

This repository provides a reference Python implementation of the explicit functional formulas developed in the paper

> *Explicit Functional Formulas for Restricted Coprime Representations of Even Integers*  
> Andrés M. Salazar

The code implements the closed functional expressions for the counting function $g(2n, p)$ (see Theorem 4.13) and includes a built-in computational verifier that directly compares the theoretical values against brute-force enumeration.

---

## Purpose

The script serves two main purposes:

1. **Exact evaluation of the functional counting formulas** derived in the paper.
2. **Independent verification of correctness**, by comparing the functional result with a direct computational count.

This provides a fully reproducible computational validation of the theoretical framework.

---

## Requirements

- Python **3.8+**
- No external dependencies (only standard library modules: `math`, `functools`)

---
cd coprime_count

