# GitHub Copilot Instructions for MELD

This document provides guidance for AI-assisted development on the MELD codebase.

## Project Overview

MELD (Modeling Employing Limited Data) is a Python/C++/CUDA package built on top of OpenMM for
running enhanced sampling molecular dynamics simulations guided by experimental data. The codebase
has two main components:

- **Python package** (`meld/`): system building, restraints, replica exchange, analysis
- **OpenMM plugin** (`plugin/`): C++/CUDA code implementing MELD forces and integrators

Most contributions touch only the Python package.

## Repository Layout

```
meld/
  system/
    builders/         # Force-field builders (amber, etc.)
    restraints.py     # All restraint type definitions
  runner/
    transform/
      restraints/     # Restraint transforms applied each step
  remd/               # Replica exchange logic
  vault.py            # I/O (reading/writing simulation data)
  test/               # Unit and integration tests
    test_functional/  # Tests requiring a full OpenMM system
  slow_tests/         # GPU-required tests (not run in CI)
plugin/               # C++/CUDA OpenMM plugin
docs/                 # Sphinx documentation source
devtools/ci/          # CI environment files and scripts
```

## Code Style

- **Formatting**: Use [black](https://github.com/ambv/black) with default settings. Run `black meld/`
  before every commit. Do not manually format code — just let black handle it.
- **Docstrings**: Use [Google style](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html)
  for all docstrings. All new functions, classes, and methods must have full docstrings.
- **Type annotations**: All new code must include `mypy`-compatible type annotations. Run
  `mypy meld` to check. CI enforces zero type errors.

## Testing

- Run the test suite with:
  ```
  python -m unittest discover meld.test
  ```
- Unit tests go in `meld/test/`. Integration tests requiring a full OpenMM setup go in
  `meld/test/test_functional/`. Slow or GPU-required tests go in `meld/slow_tests/`.
- All new features and bug fixes must include tests. Aim for high coverage on any new logic,
  especially in restraint math and system-building code.
- Do not remove or skip existing tests.

## Type Checking

```
mypy meld
```

All new code must pass `mypy` cleanly. Adding type annotations to existing un-annotated code
is also a welcome contribution (see issue #48).

## Documentation

Follow the [Grand Unified Theory of Documentation](https://documentation.divio.com) structure:

- **Tutorials**: step-by-step learning guides for new users
- **How-to guides**: recipes for accomplishing specific real-world tasks
- **Explanations**: conceptual background on MELD's approach
- **Reference**: generated from source docstrings via Sphinx

Build the docs locally with:
```
cd docs && make html
```

## Adding New Features

### New Restraint Type
1. Define the restraint class in `meld/system/restraints.py`
2. Add the corresponding transform in `meld/runner/transform/restraints/`
3. Add the force implementation to the C++/CUDA plugin (`plugin/`)
4. Write unit tests in `meld/test/`
5. Document with Google-style docstrings and type annotations

### New Force-Field Builder
- Follow the pattern in `meld/system/builders/` (e.g., the amber builder)
- Implement the same interface as existing builders
- Add tests and docstrings

## CI Pipeline

CI runs on GitHub Actions (`.github/workflows/CI.yml`) and tests:

- Python 3.10 and 3.11 with CUDA 11.8 on Ubuntu
- `python -m unittest discover meld.test`
- `mypy meld`
- Docs build (separate job, deploys to S3 on merge to `master`)

All PRs target the `master` branch. CI must pass before merge.

## Pre-Commit Checklist

Before opening a pull request:

1. `black meld/` — auto-format all changed Python files
2. `python -m unittest discover meld.test` — all tests pass
3. `mypy meld` — no type errors
4. `cd docs && make html` — docs build cleanly (if docs were changed)
5. Update `CHANGELOG.md` with a summary of changes

## Dependency Management

- Prefer existing libraries. Only add new dependencies if absolutely necessary.
- Core dependencies are listed in `setup.py` (`install_requires`).
- Conda environment specs for CI are in `devtools/ci/gh-actions/conda-envs/`.
- Optional dependencies: `gamd-openmm`, `mpi4py` (for parallel replica exchange).

## Good First Contributions

- **Docstrings**: Many existing functions have incomplete or non-Google-style docstrings (issue #48)
- **Type annotations**: Many older functions lack `mypy`-compatible type hints
- **Tests**: Coverage is incomplete, especially for restraint types
- **Documentation**: Tutorials and how-to guides are always needed

## Releasing (Maintainers Only)

1. Commit all changes and ensure tests pass
2. Update `CHANGELOG.md`
3. Run `bumpversion patch` to increment the version
4. Push commits and tags:
   ```
   git push
   git push --tags
   ```
