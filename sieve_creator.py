import numpy as np
import math

def sieve_primes(n):
    """
    Sieve of Eratosthenes.
    Returns a boolean array is_prime[0..n].
    """
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[:2] = False

    limit = int(math.isqrt(n))
    for i in range(2, limit + 1):
        if is_prime[i]:
            is_prime[i*i : n+1 : i] = False

    return is_prime

def save_sieve(is_prime, filename="sieve.npy"):
    np.save(filename, is_prime)
    print(f"Sieve saved to: {filename}")

if __name__ == "__main__":
    N = 1_000_000
    sieve = sieve_primes(N)
    save_sieve(sieve)
