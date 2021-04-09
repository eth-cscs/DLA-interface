# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class DlaInterface(CMakePackage):
    """DLA-Interface library: Distributed Linear Algebra Interface"""

    homepage = "https://github.com/eth-cscs/DLA-Interface/wiki"
    git = "https://github.com/eth-cscs/DLA-Interface"

    maintainers = ["albestro"]

    version("develop", branch="master")

    variant("miniapps", default=False, description="Build miniapps.")
    variant("fortran", default=False, description="Build with Fortran bindings.")

    depends_on("cmake@3.14:", type="build")

    depends_on("hwloc")
    # TODO OPENMP
    depends_on("mpi")
    depends_on("scalapack")

    # TODO Backends
    # depends_on("dlaf", when="+dlaf")
    # depends_on("dplasma", when="+dplasma")
    # depends_on("elpa", when="+elpa")

    def cmake_args(self):
        spec = self.spec

        # BLAS/LAPACK/SCALAPACK
        if "^mkl" in spec:
            args = [self.define("DLAI_WITH_MKL", True)]
        else:
            args = [
                self.define("DLAI_WITH_MKL", False),
                self.define(
                    "LAPACK_LIBRARY",
                    " ".join(
                        [spec[dep].libs.ld_flags for dep in ["blas", "lapack"]]
                    ),
                    ),
                self.define("SCALAPACK_LIBRARY", spec["scalapack"].libs.ld_flags),
            ]

        args.append(self.define_from_variant("DLAI_WITH_FORTRAN", "fortran"))
        args.append(self.define_from_variant("DLAI_BUILD_MINIAPPS", "miniapps"))

        # TESTs
        args.append(self.define("DLAI_BUILD_TESTING", self.run_tests))

        return args
