import matplotlib.pyplot as plt
from coprime_count import g

n_max = 10000
primes = [5, 7, 11]

plt.figure(figsize=(10, 6))

for p in primes:
    xs = []
    ys = []
    for par in range(2, n_max + 2, 2):
        xs.append(par)
        ys.append(g(par, p))
    plt.scatter(xs, ys, s=1, alpha=0.6, label=f"$p={p}$")

plt.legend(markerscale=5)
plt.title("Coprime decomposition diagram")
plt.xlabel(r"$2n$")
plt.ylabel(r"$g(2n,p)$")
plt.grid(alpha=0.2)
plt.tight_layout()
plt.show()
