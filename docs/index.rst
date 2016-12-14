**tag**: genome annotation analysis in Python!
==============================================

**tag** is a free open-source software package for analyzing genome annotation data.
It is developed as a reusable library with a focus on ease of use.
**tag** is implemented in pure Python (no compiling required) with minimal dependencies!

.. toctree::
   :maxdepth: 1

   install
   formats
   api
   acknowledgements


What problem does **tag** solve?
--------------------------------

    | *Computational biology is 90% text formatting and ID cross-referencing!*
    |     -- discouraged graduate students everywhere

Most GFF parsers will load data into memory for you--the trivial bit--but will not group related features for you--the useful bit.
**tag** represents related features as a *feature graph* (a directed acyclic graph) which can be easily traversed and inspected.

.. code:: python

   # Calculate number of exons per gene
   for gene in gff3reader:
       exons = [subfeat for subfeat in gene if subfeat.type == 'exon']
       print('num exons:', len(exons))

See :doc:`the primer on annotation formats <formats>` for more information.


Summary
-------

The **tag** library is built around the following features:

* **parsers and writers** for reading and printing annotation data in GFF3 format (with intelligent gzip support)
* **data structures** for convenient handling of various types of GFF3 entries: annotated sequence features, directives and other metadata, embedded sequences, and comments
* **generator functions** for a variety of common and useful annotation processing tasks, which can be easily composed to create streaming pipelines
* a unified **command-line interface** for executing common processing workflows
* a stable, documented **Python API** for interactive data analysis and building custom workflows


Development
-----------

Development of the **tag** library is currently a one-man show, but I would heartily welcome contributions.
The development repository is at https://github.com/standage/tag.
Please feel free to submit comments, questions, support requests to the `GitHub issue tracker <https://github.com/standage/tag/issues>`_, or (even better) a pull request!


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
