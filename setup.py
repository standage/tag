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


setup(name='aeneas',
      version='0.1.0',
      description=('Genome annotation data analysis and management implemented'
                   ' in pure Python'),
      url='http://github.com/standage/aeneas',
      author='Daniel Standage',
      author_email='daniel.standage@gmail.com',
      license='BSD-3',
      packages=['aeneas'],
      zip_safe=True)
