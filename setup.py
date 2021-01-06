#!/usr/bin/env python
# Sham for missing PEP 517 support
import os

import setuptools

# Since the rxn_chemutils dependency requires an environment variable, the
# install_requires variable must be set here instead of setup.cfg.
setuptools.setup(
    install_requires=[
        'click >= 7.1.2',
        'numpy >= 1.19.0',
        'pandas >= 1.1.4',
        'tabulate >= 0.8.7',
        'xxhash >= 2.0.0',
        'rxn_chemutils '
        '@ git+https://{}@github.ibm.com/rxn/rxn_chemutils@latest'.format(os.environ['GHE_TOKEN']),
    ]
)
