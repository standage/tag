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


d = 'Genome annotation data analysis and management implemented in pure Python'
setup(name='tag',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description=d,
      url='http://tag.readthedocs.io',
      author='Daniel Standage',
      author_email='daniel.standage@gmail.com',
      license='BSD-3',
      packages=['tag', 'tag.cli'],
      entry_points={'console_scripts': ['tag = tag.__main__:main']},
      install_requires=['intervaltree'],
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
