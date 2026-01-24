# -*- coding: utf-8 -*-
#--------------------------------------------------------------------
# The code implements the closed functional expressions for the counting 
# function g(2n,p)(see Theorem 4.13 in the paper) and includes a built-in 
# computational verifier that directly compares the theoretical values 
# against brute-force enumeration.
#----------------------------------------------------------------------
    
import math
from functools import lru_cache

# ---------------------------------------------------------------------
# Validaciones básicas
# ---------------------------------------------------------------------

def _check_int(name, x):
    if not isinstance(x, int):
        raise TypeError(f"{name} debe ser un entero (int).")


def _check_positive_int(name, x):
    _check_int(name, x)
    if x <= 0:
        raise ValueError(f"{name} debe ser un entero positivo.")


def _check_p_prime_like(p):
    """
    El paper asume p primo >= 5.
    Blindaje mínimo: forzar p >= 5.
    """
    _check_int("p", p)
    if p < 5:
        raise ValueError("p debe ser un entero >= 5 (idealmente primo).")


# ---------------------------------------------------------------------
# Operadores y funciones básicas
# ---------------------------------------------------------------------

def delta_p(a, p):
    """Resto canónico: δ_p(a) ∈ {0,...,p-1}."""
    _check_positive_int("p", p)
    _check_int("a", a)
    return a - (a // p) * p


def H(x):
    """Heaviside: 1 si x>=0, 0 si x<0."""
    _check_int("x", x)
    return 1 if x >= 0 else 0


def D(x):
    """Discriminador: 0 si x=0, 1 si x!=0."""
    _check_int("x", x)
    return 0 if x == 0 else 1


# ---------------------------------------------------------------------
# Parámetros a(p), b(p) y selector f(n, p)
# ---------------------------------------------------------------------

@lru_cache(maxsize=None)
def M_p(p):
    """
    Para p primo >=5: δ_6(p) ∈ {1,5} y M_p(p) ∈ {1,0}.
    M=1 si δ_6(p)=1, M=0 si δ_6(p)=5.
    """
    _check_p_prime_like(p)
    r6 = delta_p(p, 6)
    if r6 not in (1, 5):
        raise ValueError(
            f"Se esperaba δ_6(p) ∈ {{1,5}} para p primo>=5, pero δ_6({p})={r6}."
        )
    return (5 - r6) // 4  # 1 si r6=1, 0 si r6=5


@lru_cache(maxsize=None)
def a_p(p):
    """Solución mínima a(p) de δ_p(6x+1)=0."""
    _check_p_prime_like(p)
    Mp = M_p(p)
    return ((p - 1) // 6) * Mp + ((5 * p - 1) // 6) * (1 - Mp)


@lru_cache(maxsize=None)
def b_p(p):
    """Solución mínima b(p) de δ_p(6x+5)=0."""
    _check_p_prime_like(p)
    return p - 1 - a_p(p)


def f(n, p):
    """
    Selector entre a(p) y b(p) según δ_3(n):
    - si δ_3(n)=1 -> elige a(p)
    - si δ_3(n)=2 -> elige b(p)
    """
    _check_int("n", n)
    _check_p_prime_like(p)
    r3 = delta_p(n, 3)
    if r3 == 0:
        raise ValueError("f(n,p) solo está definido para δ_3(n)∈{1,2}.")
    # r2 = δ_2(r3): 1 si r3=1, 0 si r3=2
    r2 = delta_p(r3, 2)
    return a_p(p) * r2 + b_p(p) * (1 - r2)


# ---------------------------------------------------------------------
# Funciones κ, \barκ, τ, ν, \barν, λ y \barλ
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
    return (1 - delta_p(r, 2)) * D(2 * a - r)


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
    strict=True fuerza integridad del cociente.
    """
    _check_int("r", r)
    _check_int("a", a)
    _check_int("b", b)
    num = H(r - a) + H(r - b)
    den = 2 - D(a + b - r)
    if strict and (num % den != 0):
        raise ArithmeticError(
            f"Cociente no entero en λ: num={num}, den={den}, (r,a,b)=({r},{a},{b})"
        )
    return (r + 1) - (num // den)


def lam_bar(r, a, b, p, *, strict=True):
    """
    \barλ(r,a,b,p) = p-r-1 - (H(a-r-1)+H(b-r-1)) / (2 - D(a+b-p-r))
    strict=True fuerza integridad del cociente.
    """
    _check_int("r", r)
    _check_int("a", a)
    _check_int("b", b)
    _check_p_prime_like(p)
    num = H(a - r - 1) + H(b - r - 1)
    den = 2 - D(a + b - p - r)
    if strict and (num % den != 0):
        raise ArithmeticError(
            f"Cociente no entero en \\barλ: num={num}, den={den}, (r,a,b,p)=({r},{a},{b},{p})"
        )
    return (p - r - 1) - (num // den)


# ---------------------------------------------------------------------
# Funciones m(n), η, \barη
# ---------------------------------------------------------------------

def h_poly(x):
    """h(0)=3, h(1)=1, h(2)=5."""
    _check_int("x", x)
    return 3 * x * x - 5 * x + 3


def m(n):
    """m(n) := (n - h(δ_3(n))) / 3."""
    _check_int("n", n)
    r3 = delta_p(n, 3)
    return (n - h_poly(r3)) // 3


def m0(n):
    """Caso δ_3(n)=0."""
    _check_int("n", n)
    return (n - 3) // 3


def eta(mn, p):
    """η(m,p) := (ω+1, floor(ω/2)+1), con ω = floor(m/p)."""
    _check_int("mn", mn)
    _check_p_prime_like(p)
    wp = (mn - delta_p(mn, p)) // p
    return (wp + 1, (wp // 2) + 1)


def eta_bar(mn, p):
    """\\barη(m,p) := (ω, floor((ω-1)/2)+1)."""
    _check_int("mn", mn)
    _check_p_prime_like(p)
    wp = (mn - delta_p(mn, p)) // p
    return (wp, ((wp - 1) // 2) + 1)


# ---------------------------------------------------------------------
# Conteo funcional Q(n,p) y Q0(n,p)
# ---------------------------------------------------------------------

def Q(n, p, *, strict=True):
    """Q(n,p) aplica solo cuando δ_3(n)∈{1,2}."""
    _check_int("n", n)
    _check_p_prime_like(p)
    if delta_p(n, 3) == 0:
        raise ValueError("Q(n,p) solo aplica cuando δ_3(n)∈{1,2}. Use Q0 o Q_total.")

    mn = m(n)
    r = delta_p(mn, p)
    alpha = f(n, p)

    X1 = eta(mn, p)
    X2 = eta_bar(mn, p)
    Y1 = nu(r, alpha)
    Y2 = nu_bar(r, alpha, p)

    return X1[0] * Y1[0] + X1[1] * Y1[1] + X2[0] * Y2[0] + X2[1] * Y2[1]


def Q0(n, p, *, strict=True):
    """Q0(n,p) aplica solo cuando δ_3(n)=0."""
    _check_int("n", n)
    _check_p_prime_like(p)
    if delta_p(n, 3) != 0:
        raise ValueError("Q0(n,p) solo aplica cuando δ_3(n)=0. Use Q o Q_total.")

    mn0 = m0(n)
    r = delta_p(mn0, p)
    ap = a_p(p)
    bp = b_p(p)

    w0 = (mn0 - delta_p(mn0, p)) // p  # omega0(n,p), calculado una sola vez
    return (w0 + 1) * lam(r, ap, bp, strict=strict) + w0 * lam_bar(r, ap, bp, p, strict=strict)


def Q_total(n, p, *, strict=True):
    """Teorema total."""
    _check_int("n", n)
    _check_p_prime_like(p)
    return Q0(n, p, strict=strict) if delta_p(n, 3) == 0 else Q(n, p, strict=strict)


# ---------------------------------------------------------------------
# Conteo computacional explícito (h<=k)
# ---------------------------------------------------------------------

def contar_parejas(entrada_2n, p):
    """Conteo directo: número de pares (h,k) con h+k=2n, h<=k, gcd(h,6p)=gcd(k,6p)=1."""
    _check_int("entrada_2n", entrada_2n)
    _check_p_prime_like(p)
    if entrada_2n < 2 or entrada_2n % 2 != 0:
        raise ValueError("entrada_2n debe ser par y >= 2.")

    mod = 6 * p
    total = 0
    for h in range(1, entrada_2n // 2 + 1):
        k = entrada_2n - h
        if math.gcd(h, mod) == 1 and math.gcd(k, mod) == 1:
            total += 1
    return total


# ---------------------------------------------------------------------
# Verificador puntual
# ---------------------------------------------------------------------

def verificar(entrada_2n, p=5, *, verbose=True, devolver_detalles=True, strict=True):
    """
    Compara teoría vs cómputo directo.

    - strict=True: activa blindajes aritméticos.
    - devolver_detalles=False: retorna (teoria, computo)
    - devolver_detalles=True: retorna dict con intermedios
    """
    _check_int("entrada_2n", entrada_2n)
    _check_p_prime_like(p)
    if entrada_2n < 2 or entrada_2n % 2 != 0:
        raise ValueError("entrada_2n debe ser par y >= 2.")

    n = entrada_2n // 2
    teoria = Q_total(n, p, strict=strict)
    computo = contar_parejas(entrada_2n, p)
    ok = (teoria == computo)

    if not devolver_detalles:
        if verbose:
            print(f"2n = {entrada_2n}, p = {p}")
            print(f"Teoría:  {teoria}")
            print(f"Cómputo: {computo}")
            print("OK" if ok else "Discrepancia")
        return teoria, computo

    r3 = delta_p(n, 3)
    mn = m0(n) if r3 == 0 else m(n)
    r = delta_p(mn, p)

    out = {
        "2n": entrada_2n,
        "n": n,
        "p": p,
        "caso_delta3(n)": r3,
        "m": mn,
        "delta_p(m)": r,
        "a(p)": a_p(p),
        "b(p)": b_p(p),
        "teoria": teoria,
        "computo": computo,
        "diferencia": teoria - computo,
        "ok": ok,
        "strict": strict,
    }

    if verbose:
        print(f"2n = {entrada_2n}, p = {p}")
        print(f"caso: δ3(n) = {r3}, m = {mn}, δp(m) = {r}")
        print(f"Teoría:  {teoria}")
        print(f"Cómputo: {computo}")
        print("OK" if ok else "Discrepancia")

    return out


if __name__ == "__main__":
    print(verificar(40086, p=67, devolver_detalles=False))
