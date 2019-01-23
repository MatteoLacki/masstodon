# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='masstodon',
    packages=find_packages(),
    version='0.15',
    description='Investigate mass spectra for chemical substances, especially ETD products.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=u'Mateusz Krzysztof Łącki',
    author_email='matteo.lacki@gmail.com',
    url='https://github.com/MatteoLacki/masstodon',
    # download_url='https://github.com/MatteoLacki/masstodon/tree/GutenTag',
    keywords=[
        'Mass Spectrometry',
        'Mass spectra annotation.'
        'Analytical Chemistry',
        'ETD',
        'Electron Transfer Dissociation',
        'Fragmentation'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.5'],
    install_requires=[
        'numpy',
        'scipy',
        'lxml',
        'pyteomics>=3.4.1',
        'cffi',
        'IsoSpecPy==1.0.7',
        'networkx>=2.0',
        'matplotlib',
        'intervaltree'],
    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    # $ pip install -e ./[dev,test/] # for ZSH
    extras_require={
        'dev':  ['ipython', 'matplotlib', 'plotly'],
        'test': ['pytest']
    },
    scripts=[
        'bin/plotly_spectrum',
        'bin/masstodon'
        # 'bin/masstodon_example_call',
        # 'bin/plot_mass_spectrum',
        # 'bin/json2masstodon'
    ]
)
