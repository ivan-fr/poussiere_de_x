# 301 Smale-17 diagnostic benchmark

This is empirical only: no proof of average polynomial complexity, no external homotopy baseline.

Cases: [(2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 8), (2, 10), (3, 2), (3, 3), (3, 4), (3, 5), (4, 2), (4, 3), (5, 2), (6, 2)]
Seeds: [0, 1, 2, 3, 4, 5, 6, 7]
Base args: `--pool 4096 --epochs 32 --accept 1e-8 --tol 1e-12 --line-search 12 --probe-candidates 10 --trial-timeout 0 --keep-trials 160`
Requested roots: `min(Bezout, 20)`

| case | runs OK | roots | success runs | median sec/root | median trials/root | median failures/root | median residual | worst median residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ks(2,2) | 8/8 | 32/32 | 8/8 | 0.017530709505081177 | 7.25 | 2.375 | 3.37478383758153e-15 | 4.360393982518402e-09 |
| ks(2,3) | 8/8 | 72/72 | 8/8 | 0.03574289215935601 | 10.166666666666668 | 2.5 | 8.304602006099249e-14 | 9.257115758712225e-13 |
| ks(2,4) | 8/8 | 128/128 | 8/8 | 0.06041094660758972 | 15.1875 | 2.71875 | 8.726585506113123e-14 | 4.591667681496435e-13 |
| ks(2,5) | 8/8 | 160/160 | 8/8 | 0.03546634912490845 | 7.875 | 1.5 | 1.6549568604058202e-13 | 6.958967189061294e-13 |
| ks(2,6) | 8/8 | 160/160 | 8/8 | 0.025966376066207886 | 4.675 | 1.05 | 1.3460850028476744e-13 | 1.7304342048673278e-12 |
| ks(2,8) | 8/8 | 160/160 | 8/8 | 0.02567855715751648 | 5.6 | 2.75 | 7.805964587838033e-14 | 4.0070520727511793e-13 |
| ks(2,10) | 8/8 | 160/160 | 8/8 | 0.024967807531356814 | 5.75 | 3.125 | 1.1278361379725316e-13 | 1.920990920760983e-12 |
| ks(3,2) | 8/8 | 64/64 | 8/8 | 0.027518436312675476 | 9.125 | 2.4375 | 1.4675538099574941e-12 | 2.0876094294196747e-11 |
| ks(3,3) | 8/8 | 160/160 | 8/8 | 0.017562043666839597 | 3.425 | 0.35 | 4.0794766399925736e-14 | 5.893651107648941e-13 |
| ks(3,4) | 8/8 | 160/160 | 8/8 | 0.015448760986328126 | 2.25 | 0.175 | 3.0527870122392327e-14 | 1.5059651145692972e-12 |
| ks(3,5) | 8/8 | 160/160 | 8/8 | 0.015412348508834838 | 2.175 | 0.225 | 6.610884214362931e-14 | 6.184265406614608e-12 |
| ks(4,2) | 8/8 | 128/128 | 8/8 | 0.08251989632844925 | 23.65625 | 4.6875 | 5.5606282419317054e-14 | 1.295306235854445e-10 |
| ks(4,3) | 8/8 | 160/160 | 8/8 | 0.01354609727859497 | 2.15 | 0.225 | 2.353284574174441e-14 | 6.774726086115216e-14 |
| ks(5,2) | 8/8 | 160/160 | 8/8 | 0.01498739719390869 | 3.45 | 0.65 | 2.2316133391303236e-14 | 1.3652520467976675e-13 |
| ks(6,2) | 8/8 | 160/160 | 8/8 | 0.01219414472579956 | 2.2 | 0.25 | 3.281061905024705e-14 | 2.558089972720477e-13 |

## Global aggregate

- Runs OK: 120/120
- Roots: 2024/2024
- Successful runs: 120/120
- Median seconds/root: 0.02127094864845276
- Median trials/root: 4.675
- Median failures/root: 1.25
- Median of run median residuals: 5.895612272751972e-14
- Median evals/root: 1262.775
- Median slopes/root: 343.65

## Smale-17 candidacy notes

- Positive empirical sign: full coverage on a scaling Kostlan suite supports local robustness as a research prototype.
- Negative empirical sign: this benchmark only asks for a capped number of roots, not certified approximate zeros for arbitrary input nor all roots.
- Missing theory: no average polynomial complexity proof, no alpha/gamma/mu condition analysis, no probability bound for reaching a local cubic basin, no external homotopy comparison.