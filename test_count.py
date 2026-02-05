from coprime_count import g, count_bruteforce
import time

# WARNING: brute force on large 2n can be extremely slow.
n_test = 10000
prime = 13

# Force n_test to be even for semantic consistency
n_test = (n_test // 2) * 2

print(f"Testing g(2n, {prime}) up to 2n = {n_test}")
print("=" * 60)

start_time = time.time()
total_tests = n_test // 2

for idx, par in enumerate(range(2, n_test + 2, 2), start=1):
    gt = g(par, prime)
    bf = count_bruteforce(par, prime)
    
    if gt != bf:
        print(f"\n✗ MISMATCH at 2n={par}, p={prime}: g={gt}, brute={bf}, diff={gt-bf}")
        break
    
    # Progress indicator every 10%
    if idx % (total_tests // 10) == 0 or idx == total_tests:
        elapsed = time.time() - start_time
        progress = (idx / total_tests) * 100
        print(f"Progress: {progress:5.1f}% ({idx:5d}/{total_tests}) | "
              f"Elapsed: {elapsed:6.1f}s | "
              f"Current 2n={par}")
else:
    elapsed = time.time() - start_time
    print("=" * 60)
    print(f"✓ All tests passed! ({total_tests} tests in {elapsed:.1f}s)")
    print(f"All good up to 2n={n_test}, p={prime}")






