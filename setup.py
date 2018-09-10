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

setup(
    name='masstodon',
    packages=find_packages(),
    version='0.11',
    description='Investigate mass spectra for chemical substances, especially ETD products.',
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
        'Programming Language :: Python :: 3.6'],
    install_requires=[
        'numpy',
        'scipy',
        # 'pyteomics>=3.4.1',
        # 'lxml',
        'cvxopt',
        'cffi',
        'IsoSpecPy<2.0',
        'networkx>=2.0',
        'matplotlib',
        'intervaltree'
        ],
    # scripts=[
    #     'bin/masstodon',
    #     'bin/masstodon_example_call',
    #     'bin/plot_mass_spectrum',
    #     'bin/json2masstodon'
    # ]
)
