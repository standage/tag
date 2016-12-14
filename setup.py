#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from setuptools import setup
import glob
import versioneer


setup(name='tag',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Toolkit for genome annotation data analysis',
      url='http://github.com/standage/tag',
      author='Daniel Standage',
      author_email='daniel.standage@gmail.com',
      license='BSD-3',
      packages=['tag'],
      tests_require=['pytest', 'pytest-cov', 'pep8'],
      scripts=list(glob.glob('scripts/*.py')),
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      zip_safe=True)
