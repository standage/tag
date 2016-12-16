Installing and using **tag**
============================

The easiest way to install **tag** is with the :code:`pip` command.
By default, :code:`pip` installs the most recent stable release from the Python Package Index (PyPI).

.. code::

   pip install tag

To install the latest unreleased code, install directly from the development repository at GitHub.

.. code::

   pip install git+https://github.com/standage/tag.git

.. note:: We recommend installing **tag** and its prerequisites in a `virtual environment <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_.

**tag** is implemented in pure Python and requires no compilation.
It has only a single runtime dependency (the `intervaltree library <https://pypi.python.org/pypi/intervaltree>`_) which is also pure Python.


Using **tag** interactively
---------------------------

If you want to analyze or explore your data interactively, fire up the Python interpreter by invoking the :code:`python` command in your terminal.
Please see the :doc:`API documentation <api>` for a description of the data structures and objects available to you.

.. code:: python

   >>> import tag
   >>> reader = tag.reader.GFF3Reader(infilename='/data/genomes/mybug.gff3.gz')
   >>> for entry in tag.select.features(reader, type='intron'):
   ...     if len(entry) > 100000:
   ...         print(entry.slug)
   intron@scaffold3[37992, 149255]
   intron@scaffold55[288477, 389001]
   intron@scaffold192[1057, 196433]


Using the **tag** command-line interface
----------------------------------------

The **tag** package has a command-line interface for common processing workflows.
Execute :code:`tag -h` to see a list of available commands and :code:`tag <cmd> -h` for instructions on running a particular command.
