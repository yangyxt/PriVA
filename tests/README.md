# Test source and fixed inputs

This folder contains regression-test source and small fixed input fixtures for
installer routing, annotation interfaces, cache readers and scientific functions.
`fixtures/` includes frozen reference source used by comparisons; it is not a
place for storing completed run results.

| Source group | Relationship to PriVA |
| --- | --- |
| `test_install_utils_*.py`, `test_installer_*.py` | Exercise `scripts/install_utils.sh` and its configuration/resource helpers |
| Annotation/cache tests | Exercise the corresponding installer helpers and runtime readers |
| Scientific-function tests | Exercise the referenced `scripts/acmg_*.py` and mechanism functions |

Tests are not invoked by `main_install` or normal family prioritization.
Temporary files, generated pickles, test reports, benchmarks and execution logs
must be written outside the checkout. Use the temporary/output options provided
by each driver. See the [user manual](../README.md) for installation and runtime use.
