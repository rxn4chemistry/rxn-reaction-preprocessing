#!/usr/bin/env python

import os
from setuptools import setup, find_packages


if os.path.exists("README.md"):
    long_description = open("README.md").read()
else:
    long_description = """SMILES standardizer (e.g. metal disconnector)"""

setup(
    name="data_preprocessor",
    version="0.2.0",
    author="Alessandra Toniato, Daniel Probst @ IBM",
    author_email="ato@zurich.ibm.com",
    description="Preprocessing package for RXN Reactions",
    long_description=long_description,
    license="MIT",
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    packages=find_packages(),
    scripts=["bin/rxn-preprocess"],
)
