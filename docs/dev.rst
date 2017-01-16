Development documentation
=========================

Development environment
-----------------------

.. code::

   # Fork repo on GitHub then clone
   git clone git@github.com:YourGithubUsername/tag.git
   cd tag/

   # Optional (but highly recommended): set up a virtual environment
   virtualenv -p python3 env
   echo "export PYTHONPATH=$(pwd)" >> env/bin/activate
   source env/bin/activate  # Re-run this any time you return after launching a new terminal.

   # Install Python libraries needed for testing and run the test suite
   make devenv
   make test
   echo "make style" > .git/hooks/pre-commit

   # Test the development CLI; end users will invoke `tag` instead of `./tagcli`
   ./tagcli -h

Testing
-------

Complete coverage from automated tests is a high priority. This doesn't mean the
software is completely bug free, but it forces us to consider all the branching
behavior and make sure all code is at least executed. Our testing philosophy is
"stupidity-driven development":
* Test core features extensively, but don't try to be clever. Simple tests
  are sufficient and on balance are superior to complicated tests.
* When bugs appear, write regression tests, fix the bugs, and get back to
  something more important.
* Don't invest too much time or effort in edge features, refactoring,
  optimization, etc. If there is a simple and direct path to a substantial
  improvement, go for it and then move on. Priorities: accuracy, then
  performance, then user experience.

API
---

The **tag** package is under `semantic versioning <http://semver.org/>`_. Any
changes to the command-line interface or the Python interface require a version
bump: a minor version bump for backwards-compatible additions, and a major
version bump for API-breaking changes.

The `API documentation <http://tag.readthedocs.io/en/stable/index.html>`_ is
pulled directly from docstrings in the class and module definitions. Any Python
code in the API docs is executed and evaluated via autodoc when the test suite
is invoked.

Change log
----------

The `CHANGELOG.md` file follows the conventions described at
http://keepachangelog.com/. Minor changes that affect only the package internals
need not be documented. More substantial changes, or any bug fixes or changes to
the API, should be documented.
