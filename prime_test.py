import numpy as np
import math
from pathlib import Path

SIEVE_FILE = Path(__file__).parent / "sieve.npy"

_criba = None
_primos = None  # cache de primos como array/lista

def cargar_criba():
    """Carga la criba desde disco (memory-mapping). Se carga una sola vez."""
    global _criba
    if _criba is None:
        if not SIEVE_FILE.exists():
            raise FileNotFoundError(
                f"No se encontró {SIEVE_FILE}. "
                "Genere el archivo ejecutando sieve_creator.py."
            )
        _criba = np.load(SIEVE_FILE, mmap_mode="r")
    return _criba

def cargar_primos():
    """Construye y cachea la lista de primos usando la criba."""
    global _primos
    if _primos is None:
        criba = cargar_criba()
        _primos = np.flatnonzero(criba)  # índices True
    return _primos

def iter_primos():
    """Itera primos en orden creciente usando la lista cacheada."""
    for p in cargar_primos():
        yield int(p)

def es_primo(n: int) -> bool:
    """
    Determina si n es primo.

    - Acceso directo a la criba si n <= N.
    - División por primos si n > N, siempre que sqrt(n) <= N.
    """
    criba = cargar_criba()
    N = len(criba) - 1

    if n < 2:
        return False

    if n <= N:
        return bool(criba[n])

    limite = math.isqrt(n)
    if limite > N:
        raise ValueError(
            f"No puedo clasificar {n}: se requiere criba hasta {limite}, "
            f"pero solo existe hasta {N}."
        )

    for p in iter_primos():
        if p > limite:
            break
        if n % p == 0:
            return False
    return True

if __name__ == "__main__":
    print("🔍 Test de primalidad (Ctrl+C para salir)")
    while True:
        try:
            n = int(input("n = "))
            print(f"{n} → {'primo' if es_primo(n) else 'compuesto'}\n")
        except KeyboardInterrupt:
            print("\n👋 Saliendo.")
            break
        except Exception as e:
            print("⚠️ Error:", e, "\n")
