Acknowledgements
================

The initial development of the **tag** package was influenced heavily by `GenomeTools <http://genometools.org>`_: the `software <https://github.com/genometools/genometools>`_, the `paper <http://dx.doi.org/10.1109/TCBB.2013.68>`_, and the `developers <https://github.com/genometools/genometools/graphs/contributors>`_.
Aside from being an extremely well engineered software library and fostering a welcoming development community, GenomeTools is exceptional in two key respects: 1) a focus on streaming processing; and 2) a focus on explicit grouping of related genome features into graphs, rather than entry-by-entry processing requiring the user to resolve relationships between features.
For all of these reasons I use the GenomeTools library extensively in my research: I use the command-line interface, I've written programs (and entire libraries) that use and extend the C API, and I've contributed to the core library itself.

However, it's no secret that development in C is no walk in the park.
I've spent more hours troubleshooting memory management and chasing memory leaks than I ever care to again.
The performance of bare-metal ANSI C is nigh unmatchable, but as a research scientist *my ability to implement and evaluate prototypes quickly* saves me much more time in the long run than a constant-factor speedup in my research code's performance.

So I'd like to acknowledge the GenomeTools community for blazing the trail and (especially Gordon Gremme and Sascha Steinbiss) for their tireless support.
I'd also like to thank the Python community for their support of a wide variety of tools that made implementing, testing, documenting, and distributing the **tag** package a pleasure.
