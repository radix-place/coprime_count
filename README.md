# coprime_count.py

**Exact functional formulas for restricted coprime representations of even integers**

This repository provides a reference Python implementation of the explicit functional formulas developed in the paper

> *Explicit Functional Formulas for Restricted Coprime Representations of Even Integers*  
> Andrés M. Salazar

The code implements the functional framework introduced in the article, including:
- the canonical remainder operator,
- the functional Euclidean algorithm,
- the explicit minimal solutions of functional equations,
- and the closed-form counting formulas for coprime representations of even integers.

The implementation is designed for clarity, mathematical transparency, and direct verification of the theoretical results.

---

## Mathematical background

For a fixed prime $p \geq 5$, the script evaluates the exact number of pairs
$$
(h,k) \in \mathbb{Z}_{>0}^2, \qquad h + k = 2n, \qquad \gcd(h,6p)=\gcd(k,6p)=1,
$$
using the explicit functional formulas derived in the paper.

The theoretical value is compared, if desired, with a direct computational enumeration, providing a complete validation of the formulas.

---

## Requirements

- Python **3.8+**
- No external dependencies (only standard library modules: `math`, `functools`)

---

## Basic usage

### Direct execution

```bash
python coprime_count.py
