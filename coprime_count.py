# -*- coding: utf-8 -*-

import math
from functools import lru_cache
from prime_test import es_primo

# ---------------------------------------------------------------------
# Basic validations
# ---------------------------------------------------------------------

def _check_int(name, x):
    # NOTE: bool is a subclass of int in Python. If you want to reject bools, uncomment:
    # if isinstance(x, bool) or not isinstance(x, int):
    if not isinstance(x, int):
        raise TypeError(f"{name} must be an integer (int).")


def _check_positive_int(name, x):
    _check_int(name, x)
    if x <= 0:
        raise ValueError(f"{name} must be a positive integer.")


def _check_p_prime_like(p):
    """
    The paper assumes p is a prime >= 5.
    Exact validation using prime_test.py.
    """
    _check_int("p", p)
    if p < 5:
        raise ValueError("p must be a prime integer >= 5.")
    if not es_primo(p):
        raise ValueError(f"p must be prime >= 5. Received p={p} (not prime).")


# ---------------------------------------------------------------------
# Operators and basic functions
# ---------------------------------------------------------------------

def delta_p(a, p):
    """Canonical remainder: δ_p(a) ∈ {0,...,p-1}."""
    _check_positive_int("p", p)
    _check_int("a", a)
    return a % p  # faster and clearer than a - (a//p)*p


def H(x):
    """Heaviside: 1 if x>=0, 0 if x<0."""
    _check_int("x", x)
    return 1 if x >= 0 else 0


def D(x):
    """Discriminator: 0 if x=0, 1 if x!=0."""
    _check_int("x", x)
    return 0 if x == 0 else 1


# ---------------------------------------------------------------------
# Parameters a(p), b(p) and selector f(n, p)
# ---------------------------------------------------------------------

@lru_cache(maxsize=None)
def M_p(p):
    """
    For prime p >= 5: δ_6(p) ∈ {1,5} and M_p(p) ∈ {1,0}.
    M=1 if δ_6(p)=1, M=0 if δ_6(p)=5.
    """
    _check_p_prime_like(p)
    r6 = p % 6
    if r6 not in (1, 5):
        raise ValueError(
            f"Expected δ_6(p) ∈ {{1,5}} for prime p>=5, but δ_6({p})={r6}."
        )
    return (5 - r6) // 4  # 1 if r6=1, 0 if r6=5


@lru_cache(maxsize=None)
def a_p(p):
    """Minimal solution a(p) of δ_p(6x+1)=0."""
    _check_p_prime_like(p)
    Mp = M_p(p)
    return ((p - 1) // 6) * Mp + ((5 * p - 1) // 6) * (1 - Mp)


@lru_cache(maxsize=None)
def b_p(p):
    """Minimal solution b(p) of δ_p(6x+5)=0."""
    _check_p_prime_like(p)
    return p - 1 - a_p(p)


def f(n, p):
    """
    Selector between a(p) and b(p) based on δ_3(n):
    - if δ_3(n)=1 -> choose a(p)
    - if δ_3(n)=2 -> choose b(p)
    """
    _check_int("n", n)
    _check_p_prime_like(p)
    r3 = n % 3
    if r3 == 0:
        raise ValueError("f(n,p) is only defined for δ_3(n)∈{1,2}.")
    r2 = r3 % 2  # 1 if r3=1, 0 if r3=2
    return a_p(p) * r2 + b_p(p) * (1 - r2)


# -------------------------
# FAST INTERNALS (assume p already validated prime-like)
# -------------------------

def _f_fast_given_r3(r3, p):
    # r3 ∈ {1,2} already ensured by caller
    r2 = r3 & 1  # same as r3 % 2 for {1,2}
    ap = a_p(p)  # cached
    bp = b_p(p)  # cached
    return ap * r2 + bp * (1 - r2)


# ---------------------------------------------------------------------
# Functions κ, \barκ, τ, ν, \barν, λ and \barλ
# ---------------------------------------------------------------------

def kappa(r, a):
    _check_int("r", r)
    _check_int("a", a)
    return (r // 2) + 1 - H(r - a)


def kappa_bar(r, a, p):
    _check_int("r", r)
    _check_int("a", a)
    _check_p_prime_like(p)
    return ((p - r) // 2) - H(a - r - 1)


def tau(r, a):
    _check_int("r", r)
    _check_int("a", a)
    return (1 - (r % 2)) * D(2 * a - r)


def nu(r, a):
    t = tau(r, a)
    return (kappa(r, a) - t, t)


def nu_bar(r, a, p):
    _check_int("r", r)
    _check_int("a", a)
    _check_p_prime_like(p)
    t = tau(p - r, a - r)
    return (kappa_bar(r, a, p) - t, t)


def lam(r, a, b, *, strict=True):
    """
    λ(r,a,b) = r+1 - (H(r-a)+H(r-b)) / (2 - D(a+b-r))
    strict=True enforces integrality of the quotient.
    """
    _check_int("r", r)
    _check_int("a", a)
    _check_int("b", b)
    num = H(r - a) + H(r - b)
    den = 2 - D(a + b - r)
    if strict and (num % den != 0):
        raise ArithmeticError(
            f"Non-integer quotient in λ: num={num}, den={den}, (r,a,b)=({r},{a},{b})"
        )
    return (r + 1) - (num // den)


def lam_bar(r, a, b, p, *, strict=True):
    """
    \barλ(r,a,b,p) = p-r-1 - (H(a-r-1)+H(b-r-1)) / (2 - D(a+b-p-r))
    strict=True enforces integrality of the quotient.
    """
    _check_int("r", r)
    _check_int("a", a)
    _check_int("b", b)
    _check_p_prime_like(p)
    num = H(a - r - 1) + H(b - r - 1)
    den = 2 - D(a + b - p - r)
    if strict and (num % den != 0):
        raise ArithmeticError(
            f"Non-integer quotient in \\barλ: num={num}, den={den}, (r,a,b,p)=({r},{a},{b},{p})"
        )
    return (p - r - 1) - (num // den)


# FAST internal versions (assume p validated)
def _kappa_bar_fast(r, a, p):
    return ((p - r) // 2) - (1 if (a - r - 1) >= 0 else 0)


def _nu_bar_fast(r, a, p):
    # uses tau and kappa_bar internal fast
    t = (1 - ((p - r) & 1)) * (0 if (2 * (a - r) - (p - r)) == 0 else 1)
    return (_kappa_bar_fast(r, a, p) - t, t)


# ---------------------------------------------------------------------
# Functions m(n), η, \barη
# ---------------------------------------------------------------------

def h_poly(x):
    """h(0)=3, h(1)=1, h(2)=5."""
    _check_int("x", x)
    return 3 * x * x - 5 * x + 3


def m(n):
    """m(n) := (n - h(δ_3(n))) / 3."""
    _check_int("n", n)
    r3 = n % 3
    return (n - h_poly(r3)) // 3


def m0(n):
    """Case δ_3(n)=0."""
    _check_int("n", n)
    return (n - 3) // 3


def eta(mn, p):
    """η(m,p) := (ω+1, floor(ω/2)+1), with ω = floor(m/p)."""
    _check_int("mn", mn)
    _check_p_prime_like(p)
    wp = mn // p
    return (wp + 1, (wp // 2) + 1)


def eta_bar(mn, p):
    """\\barη(m,p) := (ω, floor((ω-1)/2)+1)."""
    _check_int("mn", mn)
    _check_p_prime_like(p)
    wp = mn // p
    return (wp, ((wp - 1) // 2) + 1)


# FAST internal versions (assume p validated)
def _eta_fast(mn, p):
    wp = mn // p
    return (wp + 1, (wp // 2) + 1)


def _eta_bar_fast(mn, p):
    wp = mn // p
    return (wp, ((wp - 1) // 2) + 1)


# ---------------------------------------------------------------------
# Functional count Q(n,p) and Q0(n,p)
# ---------------------------------------------------------------------

def Q(n, p, *, strict=True):
    """
    Q(n,p) applies only when δ_3(n)∈{1,2}.
    NOTE: 'strict' is kept for API compatibility; in this branch it does not affect arithmetic
    in the current formula (unlike Q0 where it guards λ and \barλ).
    """
    _check_int("n", n)
    _check_p_prime_like(p)
    if (n % 3) == 0:
        raise ValueError("Q(n,p) applies only when δ_3(n)∈{1,2}. Use Q0 or Q_total.")

    # Use the fast core (avoids repeating prime checks deep inside)
    return _Q_core(n, p, strict=strict)


def Q0(n, p, *, strict=True):
    """Q0(n,p) applies only when δ_3(n)=0."""
    _check_int("n", n)
    _check_p_prime_like(p)
    if (n % 3) != 0:
        raise ValueError("Q0(n,p) applies only when δ_3(n)=0. Use Q or Q_total.")

    return _Q0_core(n, p, strict=strict)


def Q_total(n, p, *, strict=True):
    """Total theorem."""
    _check_int("n", n)
    _check_p_prime_like(p)
    return _Q_total_core(n, p, strict=strict)


# -------------------------
# CORE (assumes p already validated prime-like)
# -------------------------

def _Q_core(n, p, *, strict=True):
    # assumes n int, p prime-like, and n%3 != 0 checked by caller OR routed by _Q_total_core
    r3 = n % 3
    if r3 == 0:
        raise ValueError("Internal _Q_core called with δ_3(n)=0. Use _Q0_core.")
    mn = m(n)  # m() uses n%3 internally; still cheap and clear
    r = mn % p
    alpha = _f_fast_given_r3(r3, p)

    X1 = _eta_fast(mn, p)
    X2 = _eta_bar_fast(mn, p)
    Y1 = nu(r, alpha)          # nu uses tau/kappa; no p-check
    Y2 = _nu_bar_fast(r, alpha, p)

    return X1[0] * Y1[0] + X1[1] * Y1[1] + X2[0] * Y2[0] + X2[1] * Y2[1]


def _Q0_core(n, p, *, strict=True):
    # assumes n int, p prime-like, and n%3 == 0 checked by caller OR routed by _Q_total_core
    mn0 = m0(n)
    r = mn0 % p
    ap = a_p(p)  # cached
    bp = b_p(p)  # cached

    w0 = mn0 // p
    return (w0 + 1) * lam(r, ap, bp, strict=strict) + w0 * lam_bar(r, ap, bp, p, strict=strict)


def _Q_total_core(n, p, *, strict=True):
    r3 = n % 3
    return _Q0_core(n, p, strict=strict) if r3 == 0 else _Q_core(n, p, strict=strict)


def g(entrada_2n, p, *, strict=True):
    """
    Theorem function:
        g(2n,p) := Q_total(n,p).
    """
    _check_int("entrada_2n", entrada_2n)
    _check_p_prime_like(p)
    if entrada_2n < 2 or (entrada_2n % 2) != 0:
        raise ValueError("2n must be even and >= 2.")

    return _Q_total_core(entrada_2n // 2, p, strict=strict)


# ---------------------------------------------------------------------
# Explicit brute-force count (h<=k)  -- MUST remain strict oracle
# ---------------------------------------------------------------------

def count_bruteforce(entrada_2n, p):
    """Direct count: number of pairs (h,k) with h+k=2n, h<=k, gcd(h,6p)=gcd(k,6p)=1."""
    _check_int("entrada_2n", entrada_2n)
    _check_p_prime_like(p)
    if entrada_2n < 2 or (entrada_2n % 2) != 0:
        raise ValueError("entrada_2n must be even and >= 2.")

    mod = 6 * p
    total = 0

    # micro-optimization: bind locals (no semantic change)
    gcd = math.gcd
    half = entrada_2n // 2

    for h in range(1, half + 1):
        k = entrada_2n - h
        if gcd(h, mod) == 1 and gcd(k, mod) == 1:
            total += 1
    return total


# ---------------------------------------------------------------------
# Pointwise verifier
# ---------------------------------------------------------------------

def check_theorem(entrada_2n, p=5, *, verbose=True, devolver_detalles=True, strict=True):
    """
    Compares theorem vs brute-force computation.

    - strict=True: enables arithmetic guards (relevant in Q0 via λ and \barλ)
    - devolver_detalles=False: returns (theorem, brute_force)
    - devolver_detalles=True: returns a dict with intermediates
    """
    _check_int("entrada_2n", entrada_2n)
    _check_p_prime_like(p)
    if entrada_2n < 2 or (entrada_2n % 2) != 0:
        raise ValueError("2n must be even and >= 2.")

    n = entrada_2n // 2

    # IMPORTANT: use core to avoid repeated primality checks down the call chain
    teoria = _Q_total_core(n, p, strict=strict)
    computo = count_bruteforce(entrada_2n, p)
    ok = (teoria == computo)

    if not devolver_detalles:
        if verbose:
            print(f"2n = {entrada_2n},  p = {p}")
            print(f"Theorem:     {teoria}")
            print(f"Brute force: {computo}")
            print("OK" if ok else "Mismatch")
        return teoria, computo

    r3 = n % 3
    mn = m0(n) if r3 == 0 else m(n)
    r = mn % p

    out = {
        "2n": entrada_2n,
        "n": n,
        "p": p,
        "case_delta3(n)": r3,
        "m": mn,
        "delta_p(m)": r,
        "a(p)": a_p(p),
        "b(p)": b_p(p),
        "theorem": teoria,
        "brute_force": computo,
        "difference": teoria - computo,
        "ok": ok,
        "strict": strict,
    }

    if verbose:
        print(f"2n = {entrada_2n},  p = {p}")
        print(f"case: δ3(n) = {r3},  m = {mn},  δp(m) = {r}")
        print(f"Theorem:     {teoria}")
        print(f"Brute force: {computo}")
        print("OK" if ok else "Mismatch")

    return out


if __name__ == "__main__":
    # WARNING: brute force on very large 2n can be extremely slow.
    check_theorem(50, p=7, devolver_detalles=False)
