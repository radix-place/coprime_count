from coprime_count import g, count_bruteforce

# WARNING: brute force on large 2n can be extremely slow.
n_test = 10000
prime = 101

# Force n_test to be even for semantic consistency
n_test = (n_test // 2) * 2

for par in range(2, n_test + 2, 2):
    gt = g(par, prime)
    bf = count_bruteforce(par, prime)
    if gt != bf:
        print(f"Wrong 2n={par}, p={prime}: g={gt}, brute={bf}, diff={gt-bf}")
        break
else:
    print(f"All good up to 2n={n_test}, p={prime}")






