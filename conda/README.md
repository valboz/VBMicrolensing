# Conda Packaging Notes (VBMicrolensing)

This folder contains the local conda recipe and maintainer notes used to support
the `conda-forge` package for VBMicrolensing.

Files:
- `conda/meta.yaml`: local recipe copy used for smoke testing in this repo
- `conda/MAINTAINING_PACKAGES.md`: concise maintainer workflow for conda-forge updates

Important:
- `conda/meta.yaml` in this repo is a maintainer/testing convenience.
- The actual published conda-forge package is built from the `vbmicrolensing-feedstock`
  repository after the package is accepted into conda-forge.
- The `test_conda_recipe` workflow validates that the recipe can build from the
  latest PyPI sdist, but it **does not publish to conda-forge.**
