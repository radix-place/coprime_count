import numpy as np
import math
from pathlib import Path

SIEVE_FILE = Path(__file__).parent / "sieve.npy"

_sieve = None
_primes = None  # cached primes as numpy array


def load_sieve():
    """Load the sieve from disk (memory-mapped). Loaded only once."""
    global _sieve
    if _sieve is None:
        if not SIEVE_FILE.exists():
            raise FileNotFoundError(
                f"{SIEVE_FILE} was not found. "
                "Generate it by running sieve_creator.py."
            )
        _sieve = np.load(SIEVE_FILE, mmap_mode="r")
    return _sieve


def load_primes():
    """Build and cache the list of primes using the sieve."""
    global _primes
    if _primes is None:
        sieve = load_sieve()
        _primes = np.flatnonzero(sieve)  # indices where sieve is True
    return _primes


def iter_primes():
    """Iterate primes in increasing order using the cached list."""
    for p in load_primes():
        yield int(p)


def is_prime(n: int) -> bool:
    """
    Determine whether n is prime.

    - Direct sieve lookup if n <= N.
    - Trial division by primes if n > N, provided sqrt(n) <= N.
    """
    sieve = load_sieve()
    N = len(sieve) - 1

    if n < 2:
        return False

    if n <= N:
        return bool(sieve[n])

    limit = math.isqrt(n)
    if limit > N:
        raise ValueError(
            f"Cannot classify {n}: sieve up to {limit} is required, "
            f"but the current sieve only goes up to {N}."
        )

    for p in iter_primes():
        if p > limit:
            break
        if n % p == 0:
            return False
    return True


if __name__ == "__main__":
    print("Primality test (Ctrl+C to exit)")
    while True:
        try:
            n = int(input("n = "))
            print(f"{n} → {'prime' if is_prime(n) else 'composite'}\n")
        except KeyboardInterrupt:
            print("\n Exiting.")
            break
        except Exception as e:
            print("Error:", e, "\n")
