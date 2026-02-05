import matplotlib.pyplot as plt
from coprime_count import g

n_max = 10000
primes = [5, 7, 11]
save_figure = False  # 

print(f"Generating coprime decomposition diagram for 2n ≤ {n_max}")
print(f"Primes: {primes}")
print("Computing values...")

plt.figure(figsize=(10, 6))

for p in primes:
    print(f"  Computing for p={p}...", end=" ")
    xs = []
    ys = []
    for par in range(2, n_max + 2, 2):
        xs.append(par)
        ys.append(g(par, p))
    plt.scatter(xs, ys, s=1, alpha=0.6, label=f"$p={p}$")
    print("✓")

plt.legend(markerscale=5)
plt.title("Coprime decomposition diagram")
plt.xlabel(r"$2n$")
plt.ylabel(r"$g(2n,p)$")
plt.grid(alpha=0.2)
plt.tight_layout()

# Guardar si save_figure = True
if save_figure:
    filename = f"coprime_diagram_n{n_max}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved as: {filename}")

plt.show()