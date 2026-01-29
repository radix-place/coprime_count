# -*- coding: utf-8 -*-

import math
from functools import lru_cache
from prime_test import is_prime


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
# Basic operators
# ---------------------------------------------------------------------

def H(x):
    return 1 if x >= 0 else 0


def D(x):
    return 0 if x == 0 else 1


# ---------------------------------------------------------------------
# Parameters a(p), b(p)
# ---------------------------------------------------------------------

@lru_cache(maxsize=None)
def M_p(p):
    r6 = p % 6
    # For primes p >= 5 we must have p ≡ 1 or 5 (mod 6).
    assert r6 in (1, 5), f"Unexpected p mod 6 = {r6}. Expected 1 or 5 for prime p>=5."
    return (5 - r6) // 4   # 1 if p ≡ 1 mod 6, 0 if p ≡ 5 mod 6


@lru_cache(maxsize=None)
def a_p(p):
    Mp = M_p(p)
    return ((p - 1) // 6) * Mp + ((5 * p - 1) // 6) * (1 - Mp)


@lru_cache(maxsize=None)
def b_p(p):
    return p - 1 - a_p(p)


# ---------------------------------------------------------------------
# κ, κ̄, τ, ν, ν̄, λ, λ̄
# ---------------------------------------------------------------------

def kappa(r, a):
    return (r // 2) + 1 - H(r - a)


def kappa_bar(r, a, p):
    return ((p - r) // 2) - H(a - r - 1)


def tau(r, a):
    return (1 - (r & 1)) * D(2 * a - r)


def nu(r, a):
    t = tau(r, a)
    return (kappa(r, a) - t, t)


def nu_bar(r, a, p):
    t = tau(p - r, a - r)
    return (kappa_bar(r, a, p) - t, t)


def lam(r, a, b):
    return (r + 2) - H(r - a) - H(r - b) - D(a + b - r)


def lam_bar(r, a, b, p):
    return p - r - H(a - r - 1) - H(b - r - 1) - D(a + b - p - r)


# ---------------------------------------------------------------------
# Auxiliary arithmetic functions
# ---------------------------------------------------------------------

def h_poly(x):
    return 3 * x * x - 5 * x + 3


def m(n):
    return (n - h_poly(n % 3)) // 3


def m0(n):
    return (n - 3) // 3


def eta(mn, p):
    w = mn // p
    return (w + 1, (w // 2) + 1)


def eta_bar(mn, p):
    w = mn // p
    return (w, ((w - 1) // 2) + 1)


# Fast internal versions 

def _eta_fast(mn, p):
    w = mn // p
    return (w + 1, (w // 2) + 1)


def _eta_bar_fast(mn, p):
    w = mn // p
    return (w, ((w - 1) // 2) + 1)


def _nu_bar_fast(r, a, p):
    # t = tau(p-r, a-r)
    t = (1 - ((p - r) & 1)) * (0 if (2 * (a - r) - (p - r)) == 0 else 1)
    # kappa_bar(r, a, p) = ((p-r)//2) - H(a-r-1)
    kb = ((p - r) // 2) - (1 if (a - r - 1) >= 0 else 0)
    return (kb - t, t)


# ---------------------------------------------------------------------
# Core functional counts
# ---------------------------------------------------------------------

def _Q_core(n, p):
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
    mn0 = m0(n)
    r = mn0 % p
    ap = a_p(p)
    bp = b_p(p)
    w0 = mn0 // p

    return (w0 + 1) * lam(r, ap, bp) + w0 * lam_bar(r, ap, bp, p)


def _Q_total_core(n, p):
    return _Q0_core(n, p) if (n % 3 == 0) else _Q_core(n, p)


# ---------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------

def g(entrada_2n, p):
    """
    Theorem:
        g(2n, p) = Q_total(n, p)
    """
    _check_int("2n", entrada_2n)
    _check_p_prime_like(p)
    if entrada_2n < 2 or (entrada_2n & 1):
        raise ValueError("2n must be even and >= 2.")

    return _Q_total_core(entrada_2n // 2, p)


# ---------------------------------------------------------------------
# Brute-force oracle
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
    Direct count:
    number of pairs (h,k) such that
        h + k = 2n,
        h ≤ k,
        gcd(h, 6p) = gcd(k, 6p) = 1.
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
        print("OK" if ok else "Mismatch")

    return theorem, bruteforce, ok


# ---------------------------------------------------------------------
# Test example
# ---------------------------------------------------------------------

if __name__ == "__main__":
    check_theorem(58, p=5)
