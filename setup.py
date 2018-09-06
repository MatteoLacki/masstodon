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
    name='MassTodonPy',
    packages=find_packages(),
    version='0.3.7',
    description='Estimate the products of Electron Transfer Dissociation in \
    Mass Spectrometry for a given biological substance and \
    the chemical reaction probabilities that lead to these products.',
    author=u'Mateusz Krzysztof Łącki',
    author_email='matteo.lacki@gmail.com',
    url='https://github.com/MatteoLacki/MassTodonPy',
    download_url='https://github.com/MatteoLacki/MassTodonPy/tree/GutenTag',
    keywords=[
        'Mass Spectrometry',
        'ETD',
        'Electron Transfer Dissociation',
        'Fragmentation'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.4'],
    install_requires=[
        'cffi',
        'numpy',
        'pyteomics>=3.4.1',
        'lxml',
        'cvxopt',
        'IsoSpecPy<2.0',
        'networkx>=2.0',
        'future',
        'six',
        'bokeh'],
    scripts=[
        'bin/masstodon',
        'bin/masstodon_example_call',
        'bin/plot_mass_spectrum',
        'bin/json2masstodon']
)
