# -*- coding: utf-8 -*-
"""
Exact Counting of Restricted Coprime Representations
Author: Andrés M. Salazar
Paper: "Exact Counting of Restricted Coprime Representations of Even Integers 
        via Functional Residue Calculus"
"""

import math
import time
from functools import lru_cache
from prime_test import is_prime

__version__ = "1.0.0"
__all__ = ['g', 'count_bruteforce', 'check_theorem', 'benchmark_comparison']


# ---------------------------------------------------------------------
# Validations 
# ---------------------------------------------------------------------

def _check_int(name, x):
    if not isinstance(x, int):
        raise TypeError(f"{name} must be an integer (int).")


def _check_p_prime_like(p):
    _check_int("p", p)
    if p < 5:
        raise ValueError("p must be a prime integer >= 5.")
    if not is_prime(p):
        raise ValueError(f"p must be prime >= 5. Received p={p} (not prime).")


# ---------------------------------------------------------------------
# Basic operators (Definition 2.1 in paper)
# ---------------------------------------------------------------------

def delta_k(a, k):
    """
    Canonical remainder operator δ_k(a).
    Returns a - k⌊a/k⌋ ∈ {0, 1, ..., k-1}.
    
    Reference: Definition 2.1 (equation 3)
    """
    return a - k * (a // k)


def H(x):
    """Step function: H(x) = 1 if x ≥ 0, else 0."""
    return 1 if x >= 0 else 0


def D(x):
    """
    Complementary indicator function:
    
        D(x) = 1 - 𝟙(x),
    
    where 𝟙(x) is the indicator of x = 0 used in the paper:
    
        𝟙(x) = 1 if x = 0, else 0.
    
    Thus:
    
        D(x) = 0 if x = 0,
        D(x) = 1 if x ≠ 0.
    
    This algebraically equivalent form is used to simplify expressions
    while preserving exact counting results.
    """
    return 0 if x == 0 else 1


# ---------------------------------------------------------------------
# Parameters a(p), b(p) (Lemmas 4.2, 4.3)
# ---------------------------------------------------------------------

@lru_cache(maxsize=None)
def M_p(p):
    """
    Selector function M(p) from equation (18).
    M(p) = 1 if δ_6(p) = 1, else 0.
    """
    r6 = p % 6
    # For primes p >= 5 we must have p ≡ 1 or 5 (mod 6).
    assert r6 in (1, 5), f"Unexpected p mod 6 = {r6}. Expected 1 or 5 for prime p>=5."
    return (5 - r6) // 4   # 1 if p ≡ 1 mod 6, 0 if p ≡ 5 mod 6


@lru_cache(maxsize=None)
def a_p(p):
    """
    Minimal solution of δ_p(6x + 1) = 0.
    Reference: Lemma 4.2, equation (18).
    """
    Mp = M_p(p)
    return ((p - 1) // 6) * Mp + ((5 * p - 1) // 6) * (1 - Mp)


@lru_cache(maxsize=None)
def b_p(p):
    """
    Minimal solution of δ_p(6x + 5) = 0.
    Reference: Lemma 4.3, equation (23).
    """
    return p - 1 - a_p(p)


# ---------------------------------------------------------------------
# κ, κ̄, τ, ν, ν̄, λ, λ̄ (Section 4.1)
# ---------------------------------------------------------------------

def kappa(r, a):
    """κ(r, a) from Lemma 4.5, equation (25)."""
    return (r // 2) + 1 - H(r - a)


def kappa_bar(r, a, p):
    """κ̄(r, a, p) from Lemma 4.7, equation (28)."""
    return ((p - r) // 2) - H(a - r - 1)


def tau(r, a):
    """τ(r, a) from equation (26)."""
    return (1 - (r & 1)) * D(2 * a - r)


def nu(r, a):
    """ν(r, a) from equation (27)."""
    t = tau(r, a)
    return (kappa(r, a) - t, t)


def nu_bar(r, a, p):
    """ν̄(r, a, p) from equation (29)."""
    t = tau(p - r, a - r)
    return (kappa_bar(r, a, p) - t, t)


def lam(r, a, b):
    """λ(r, a, b) from Lemma 4.9, equation (30)."""
    return (r + 2) - H(r - a) - H(r - b) - D(a + b - r)


def lam_bar(r, a, b, p):
    """λ̄(r, a, b, p) from Lemma 4.10, equation (31)."""
    return p - r - H(a - r - 1) - H(b - r - 1) - D(a + b - p - r)


# ---------------------------------------------------------------------
# Auxiliary arithmetic functions (Section 4.2)
# ---------------------------------------------------------------------

def h_poly(x):
    """Polynomial h(x) = 3x² - 5x + 3."""
    return 3 * x * x - 5 * x + 3


def m(n):
    """
    Auxiliary function m(n) from equation (33).
    Used when δ_3(n) ∈ {1, 2}.
    """
    return (n - h_poly(n % 3)) // 3


def m0(n):
    """
    Auxiliary function m_0(n) = (n - 3)/3.
    Used when δ_3(n) = 0.
    """
    return (n - 3) // 3


def omega(mn, p):
    """ω(n, p) from equation (33)."""
    return mn // p


def eta(mn, p):
    """η(m(n), p) from Lemma 4.11."""
    w = omega(mn, p)
    return (w + 1, (w // 2) + 1)


def eta_bar(mn, p):
    """η̄(m(n), p) from Lemma 4.11."""
    w = omega(mn, p)
    return (w, ((w - 1) // 2) + 1)


# Fast internal versions (avoid function call overhead)

def _eta_fast(mn, p):
    w = mn // p
    return (w + 1, (w // 2) + 1)


def _eta_bar_fast(mn, p):
    w = mn // p
    return (w, ((w - 1) // 2) + 1)


def _nu_bar_fast(r, a, p):
    """Optimized version of nu_bar without intermediate function calls."""
    # t = tau(p-r, a-r)
    t = (1 - ((p - r) & 1)) * (0 if (2 * (a - r) - (p - r)) == 0 else 1)
    # kappa_bar(r, a, p) = ((p-r)//2) - H(a-r-1)
    kb = ((p - r) // 2) - (1 if (a - r - 1) >= 0 else 0)
    return (kb - t, t)


# ---------------------------------------------------------------------
# Core functional counts (Lemmas 4.11, 4.12)
# ---------------------------------------------------------------------

def _Q_core(n, p):
    """
    Q(n, p) for δ_3(n) ∈ {1, 2}.
    Reference: Lemma 4.11, equation (34).
    """
    r3 = n % 3
    mn = m(n)
    r = mn % p
    alpha = a_p(p) if r3 == 1 else b_p(p)

    X1 = _eta_fast(mn, p)
    X2 = _eta_bar_fast(mn, p)
    Y1 = nu(r, alpha)
    Y2 = _nu_bar_fast(r, alpha, p)

    return X1[0] * Y1[0] + X1[1] * Y1[1] + X2[0] * Y2[0] + X2[1] * Y2[1]


def _Q0_core(n, p):
    """
    Q_0(n, p) for δ_3(n) = 0.
    Reference: Lemma 4.12.
    """
    mn0 = m0(n)
    r = mn0 % p
    ap = a_p(p)
    bp = b_p(p)
    w0 = mn0 // p

    return (w0 + 1) * lam(r, ap, bp) + w0 * lam_bar(r, ap, bp, p)


def _Q_total_core(n, p):
    """
    Combined Q function.
    Reference: Theorem 4.13.
    """
    return _Q0_core(n, p) if (n % 3 == 0) else _Q_core(n, p)


# ---------------------------------------------------------------------
# Public API (Main Theorem 4.13)
# ---------------------------------------------------------------------

def g(entrada_2n, p):
    """
    Compute g(2n, p) = number of pairs (h, k) such that:
        - h + k = 2n
        - h ≤ k
        - gcd(h, 6p) = gcd(k, 6p) = 1
    
    Reference: Theorem 4.13 (main result).
    
    Parameters
    ----------
    entrada_2n : int
        Even integer 2n ≥ 2
    p : int
        Prime number p ≥ 5
    
    Returns
    -------
    int
        Number of valid coprime pairs
    
    Examples
    --------
    >>> g(20, 5)
    2
    >>> g(100, 7)
    12
    """
    _check_int("2n", entrada_2n)
    _check_p_prime_like(p)
    if entrada_2n < 2 or (entrada_2n & 1):
        raise ValueError("2n must be even and >= 2.")

    return _Q_total_core(entrada_2n // 2, p)


# ---------------------------------------------------------------------
# Brute-force oracle (Section 5.1)
# ---------------------------------------------------------------------

def _count_bruteforce_core(entrada_2n, p):
    """Core brute-force (assumes inputs already validated)."""
    mod = 6 * p
    total = 0
    half = entrada_2n // 2
    gcd = math.gcd

    for h in range(1, half + 1):
        k = entrada_2n - h
        if gcd(h, mod) == 1 and gcd(k, mod) == 1:
            total += 1

    return total


def count_bruteforce(entrada_2n, p):
    """
    Direct enumeration count (oracle for validation).
    
    Reference: Section 5.1
    
    Parameters
    ----------
    entrada_2n : int
        Even integer 2n ≥ 2
    p : int
        Prime number p ≥ 5
    
    Returns
    -------
    int
        Number of valid coprime pairs (computed by brute force)
    """
    _check_int("2n", entrada_2n)
    _check_p_prime_like(p)
    if entrada_2n < 2 or (entrada_2n & 1):
        raise ValueError("2n must be even and >= 2.")

    return _count_bruteforce_core(entrada_2n, p)


# ---------------------------------------------------------------------
# Verifier (single validation)
# ---------------------------------------------------------------------

def check_theorem(entrada_2n, p=5, *, verbose=True):
    """
    Validate Theorem 4.13 for a single case.
    
    Parameters
    ----------
    entrada_2n : int
        Even integer 2n ≥ 2
    p : int, default=5
        Prime number p ≥ 5
    verbose : bool, default=True
        Print results
    
    Returns
    -------
    tuple (int, int, bool)
        (theorem_result, bruteforce_result, match)
    """
    _check_int("2n", entrada_2n)
    _check_p_prime_like(p)
    if entrada_2n < 2 or (entrada_2n & 1):
        raise ValueError("2n must be even and >= 2.")

    n = entrada_2n // 2
    theorem = _Q_total_core(n, p)
    bruteforce = _count_bruteforce_core(entrada_2n, p)
    ok = (theorem == bruteforce)

    if verbose:
        print(f"2n = {entrada_2n}, p = {p}")
        print(f"Theorem:     {theorem}")
        print(f"Brute force: {bruteforce}")
        print("✓ OK" if ok else "✗ MISMATCH")

    return theorem, bruteforce, ok


# ---------------------------------------------------------------------
# NEW: Benchmark and timing comparison (Section 5.3)
# ---------------------------------------------------------------------

def benchmark_single(entrada_2n, p, method='functional', warmup=3, runs=100):
    """
    Benchmark a single computation.
    
    Parameters
    ----------
    entrada_2n : int
        Even integer
    p : int
        Prime number
    method : str
        'functional' or 'bruteforce'
    warmup : int
        Number of warmup runs
    runs : int
        Number of timed runs
    
    Returns
    -------
    dict
        {'mean_time': float, 'std_time': float, 'result': int}
    """
    _check_int("2n", entrada_2n)
    _check_p_prime_like(p)
    
    n = entrada_2n // 2
    
    # Select function
    if method == 'functional':
        func = lambda: _Q_total_core(n, p)
    elif method == 'bruteforce':
        func = lambda: _count_bruteforce_core(entrada_2n, p)
    else:
        raise ValueError("method must be 'functional' or 'bruteforce'")
    
    # Warmup
    for _ in range(warmup):
        result = func()
    
    # Timing
    times = []
    for _ in range(runs):
        start = time.perf_counter()
        result = func()
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    
    import statistics
    return {
        'mean_time': statistics.mean(times),
        'std_time': statistics.stdev(times) if len(times) > 1 else 0.0,
        'result': result
    }


def benchmark_comparison(max_n=10000, p=7, step=100, runs=50):
    """
    Compare functional vs brute-force methods across range.
    
    Reference: Section 5.3
    
    Parameters
    ----------
    max_n : int
        Maximum n to test
    p : int
        Prime number
    step : int
        Step size for n
    runs : int
        Number of runs per test
    
    Returns
    -------
    dict
        Results dictionary with timing data
    """
    _check_p_prime_like(p)
    
    results = {
        'n_values': [],
        'functional_times': [],
        'bruteforce_times': [],
        'speedup': []
    }
    
    print(f"Benchmarking comparison (p={p}, max_n={max_n}, runs={runs})")
    print("-" * 60)
    
    for n in range(step, max_n + 1, step):
        entrada_2n = 2 * n
        
        # Functional method
        func_stats = benchmark_single(entrada_2n, p, 'functional', runs=runs)
        
        # Brute force
        brute_stats = benchmark_single(entrada_2n, p, 'bruteforce', runs=runs)
        
        # Verify they match
        assert func_stats['result'] == brute_stats['result'], \
            f"Mismatch at 2n={entrada_2n}"
        
        speedup = brute_stats['mean_time'] / func_stats['mean_time']
        
        results['n_values'].append(n)
        results['functional_times'].append(func_stats['mean_time'])
        results['bruteforce_times'].append(brute_stats['mean_time'])
        results['speedup'].append(speedup)
        
        if n % (step * 10) == 0:
            print(f"n={n:5d}  |  Functional: {func_stats['mean_time']*1e6:6.2f} μs  |  "
                  f"Brute: {brute_stats['mean_time']*1e6:7.2f} μs  |  "
                  f"Speedup: {speedup:6.1f}x")
    
    print("-" * 60)
    print(f"Average speedup: {sum(results['speedup'])/len(results['speedup']):.1f}x")
    
    return results


# ---------------------------------------------------------------------
# NEW: Extended validation (multiple primes, large range)
# ---------------------------------------------------------------------

def validate_range(max_n, primes, verbose_interval=1000):
    """
    Validate theorem across range of n and multiple primes.
    
    Parameters
    ----------
    max_n : int
        Maximum n to test
    primes : list of int
        List of primes to test
    verbose_interval : int
        Print progress every N tests
    
    Returns
    -------
    dict
        Validation results
    """
    results = {
        'total_tests': 0,
        'passed': 0,
        'failed': 0,
        'failures': []
    }
    
    print(f"Validating g(2n, p) for n ≤ {max_n}, primes = {primes}")
    print("-" * 60)
    
    for p in primes:
        _check_p_prime_like(p)
        for n in range(1, max_n + 1):
            entrada_2n = 2 * n
            theorem, bruteforce, ok = check_theorem(entrada_2n, p, verbose=False)
            
            results['total_tests'] += 1
            if ok:
                results['passed'] += 1
            else:
                results['failed'] += 1
                results['failures'].append((entrada_2n, p, theorem, bruteforce))
            
            if results['total_tests'] % verbose_interval == 0:
                print(f"Progress: {results['total_tests']} tests, "
                      f"{results['passed']} passed, {results['failed']} failed")
    
    print("-" * 60)
    print(f"Validation complete: {results['passed']}/{results['total_tests']} passed")
    
    if results['failed'] > 0:
        print(f"\n⚠ {results['failed']} FAILURES:")
        for entrada_2n, p, th, bf in results['failures']:
            print(f"  2n={entrada_2n}, p={p}: theorem={th}, bruteforce={bf}")
    else:
        print("✓ All tests passed!")
    
    return results


# ---------------------------------------------------------------------
# NEW: Export results to CSV
# ---------------------------------------------------------------------

def export_to_csv(filename, max_n, p):
    """
    Export g(2n, p) values to CSV for external analysis.
    
    Parameters
    ----------
    filename : str
        Output CSV filename
    max_n : int
        Maximum n
    p : int
        Prime number
    """
    import csv
    
    _check_p_prime_like(p)
    
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['2n', 'n', 'p', 'g(2n,p)', 'delta_3(n)', 'm(n)', 'delta_p(m(n))'])
        
        for n in range(1, max_n + 1):
            entrada_2n = 2 * n
            g_val = g(entrada_2n, p)
            delta3 = n % 3
            mn = m(n) if delta3 != 0 else m0(n)
            delta_p = mn % p
            
            writer.writerow([entrada_2n, n, p, g_val, delta3, mn, delta_p])
    
    print(f"Exported {max_n} values to {filename}")


# ---------------------------------------------------------------------
# Test example
# ---------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("Single verification example:")
    print("=" * 60)
    check_theorem(58, p=5)
    
    print("\n" + "=" * 60)
    print("Benchmark comparison:")
    print("=" * 60)
    benchmark_comparison(max_n=1000, p=7, step=100, runs=30)
