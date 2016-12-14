tag: genome annotation analysis in Python!
==========================================

**tag** is a free open-source software package for working with genome annotation data.
It is developed as a reusable library with a focus on ease of use.
**tag** is implemented in pure Python (no compiling required) and has no runtime dependencies!

.. toctree::
   :maxdepth: 1

   install
   gff3
   api
   acknowledgements
   naming


What problem does **tag** solve?
--------------------------------

    | *Computational biology is 90% text formatting and ID cross-referencing!*
    |     -- discouraged graduate students everywhere

Most GFF parsers will load data into memory for you--the trivial bit--but will not group related features for you--the useful bit.
**tag** represents related features as a *feature graph* (a directed acyclic graph) which can be easily traversed and inspected.

See :doc:`this page <gff3>` for more information.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
