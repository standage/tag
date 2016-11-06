#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------
"""Setup configuration for aeneas"""

from setuptools import setup
import glob


setup(name='aeneas',
      version='0.1.0',
      description='Tools for genome annotation data analysis',
      url='http://github.com/standage/aeneas',
      author='Daniel Standage',
      author_email='daniel.standage@gmail.com',
      license='BSD-3',
      packages=['aeneas'],
      tests_require=['pytest', 'pytest-cov', 'pep8'],
      scripts=list(glob.glob('scripts/*.py')),
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      zip_safe=True)
