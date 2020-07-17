#!/usr/bin/env python

import os
from setuptools import setup


if os.path.exists('README.md'):
    long_description = open('README.md').read()
else:
    long_description = '''SMILES standardizer (e.g. metal disconnector)'''

setup(
    name='data_preprocessor',
    version='0.1.0',
    author='Alessandra Toniato @ IBM',
    author_email='ato@zurich.ibm.com',
    py_modules=['data_preprocessor'],
    description='Preprocessing repo for SMILEs',
    long_description=long_description,
    license='MIT',
    requirements=[
        # 'requests==2.23.0'
        'pandas>=1.0.1',
        'crc64iso'==0.0.2,
        'action_sequences @ git+https://{}@github.ibm.com/rxn/action_sequences.git@latest#egg=action_sequences'.format(
        os.environ['GHE_TOKEN']
    ],
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ]
)
