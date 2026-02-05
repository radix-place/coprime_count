# Coprime decomposition counts and diagrams

This repository contains the reference implementation accompanying the article:

> *Exact Counting of Restricted Coprime Representations of Even Integers via Functional Residue Calculus*  
> Andrés M. Salazar  
> Pontificia Universidad Javeriana Cali — Colombia  

The code provides:
- An explicit closed formula for the function $g(2n,p)$ (see Theorem 4.13 in the paper),
- A brute-force oracle for direct verification,
- Systematic computational validation,
- Visualization scripts for generating coprime decomposition diagrams.

---

## Contents

- **`coprime_count.py`**  
  Core implementation of the closed formula for `g(2n,p)`, together with a strict brute-force oracle and pointwise verification routines.

- **`prime_test.py`**  
  Efficient primality testing based on a precomputed sieve and memory-mapped access.

- **`sieve_creator.py`**  
  Script to generate the prime sieve file `sieve.npy`.  
  *(Run this script first to create `sieve.npy`, which is required by the other scripts.)*

- **`test_count.py`**  
  Exhaustive pointwise verification of the closed formula against a direct brute-force oracle for all even integers up to a prescribed bound.

- **`graph.py`**  
  Script to generate coprime decomposition diagrams (Figure 1 in the paper).

---

## Quick Start

### 1. Setup

Install dependencies:
```bash
pip install numpy matplotlib
```

Generate the prime sieve (required for all scripts):
```bash
python sieve_creator.py
```

### 2. Basic Usage

**Single verification:**
```python
from coprime_count import g, count_bruteforce, check_theorem

# Compute g(100, 7) using the closed formula
result = g(100, 7)
print(f"g(100, 7) = {result}")

# Verify against brute force
check_theorem(100, 7)
```

**Run complete validation:**
```bash
python test_count.py
```

**Generate visualization:**
```bash
python graph.py
```

### 3. Example Output
```
Testing g(2n, 13) up to 2n = 10000
============================================================
Progress:  10.0% (  500/5000) | Elapsed:    0.0s | Current 2n=1000
Progress:  20.0% ( 1000/5000) | Elapsed:    0.0s | Current 2n=2000
...
Progress: 100.0% ( 5000/5000) | Elapsed:    1.1s | Current 2n=10000
============================================================
✓ All tests passed! (5000 tests in 1.1s)
All good up to 2n=10000, p=13
```

---

## Performance

The functional algorithm runs in **O(1)** time for fixed prime `p`, compared to **O(n)** for brute force enumeration.

Benchmark results (p=7, n=1000):
- Functional method: **~1.2 μs**
- Brute force: **~64 μs**
- **Speedup: ~53x**

---

## Requirements

- Python ≥ 3.9  
- NumPy  
- Matplotlib

---

## Citation

If you use this code in your research, please cite:
```bibtex
@article{salazar2026coprime,
  title={Exact Counting of Restricted Coprime Representations of Even Integers via Functional Residue Calculus},
  author={Salazar, Andr{\'e}s M.},
  year={2026}
}
```

---

## License

This code is released under the MIT License. See `LICENSE` for details.

---

## Contact

For questions or comments, please contact:  
**Andrés M. Salazar** — andresmsalazar@javerianacali.edu.co
