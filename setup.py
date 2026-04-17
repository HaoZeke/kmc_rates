import os

from numpy.distutils.core import setup
from numpy.distutils.core import Extension

depends = ["source/graph.hpp", "source/ngt.hpp"]
extra_compile_args = [
    "-Wall", "-Wextra", "-g", "-O2",
    "-funroll-loops", "-mtune=native", "-std=c++0x",
]
extra_link_args = []
include_dirs = ["source/"]

setup(
    name="kmc_rates",
    version="0.0.1",
    description=(
        "kmc_rates NGT reference (Rowley 2014) -- installable fork "
        "used by the lode-org/amsel cross-check fixture."
    ),
    packages=["kmc_rates"],
    ext_modules=[
        Extension(
            "kmc_rates/ngt",
            ["kmc_rates/ngt.cpp"],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
            include_dirs=include_dirs,
            language="c++",
            depends=depends,
        ),
    ],
)
