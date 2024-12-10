# LinearAlgebra

This package is part of the Julia standard library (stdlib).

`LinearAlgebra.jl` provides functionality for performing linear algebra operations in Julia.

| **Build Status**                    | **Coverage**                    |
|:-------------------------------:|:-------------------------------:|
[![Build status](https://badge.buildkite.com/4032f509b3a7354c8ce7b34b98c7496d66adc4504b0101cbcf.svg?branch=master)](https://buildkite.com/julialang/linearalgebra-dot-jl) | [![][codecov-img]][codecov-url] |

[codecov-img]: https://codecov.io/gh/JuliaLang/LinearAlgebra.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaLang/LinearAlgebra.jl

## Migrating a Pull Request (PR) from the Julia repository to this repository

Since this package was split out from the main Julia repository, you might have previously made a pull request (PR) to the Julia repo. You can easily migrate such PRs to this repository using the following steps:

1. Add the Julia repository (or your fork) as a new remote repository:
   ```bash
   git remote add juliarepo https://github.com/JuliaLang/julia
   ```

2. Fetch the commits from the Julia repository:
   ```bash
   git fetch juliarepo
   ```

3. Cherry-pick the relevant commits made in the Julia repository to this repository:
   ```bash
   git cherry-pick $JULIA_COMMIT
   ```

## Using development versions of this package

This package performs some type piracy and is also included in the sysimage, which makes using a development version slightly more complex than usual.

To use a development version of this package, you can choose one of the following methods:

1. **Change the UUID in the project file and load the package:**
   This approach will produce warnings and may lead to method ambiguities between the development version and the one in the sysimage, but it can be used for basic experimentation.

2. **Build Julia with the custom `LinearAlgebra` commit:**
   Modify the commit in `stdlib/LinearAlgebra.version` and build Julia.

3. **Build a custom sysimage with the new `LinearAlgebra`:**
   - Install `PackageCompiler`.
   - Load it and, with this project active, run:
     ```julia
     create_sysimage(["LinearAlgebra"]; sysimage_path="new_sysimage.so", incremental=false, filter_stdlibs=true)
     ```
   - Start Julia with the custom sysimage:
     ```bash
     julia -Jnew_sysimage.so
     ```
