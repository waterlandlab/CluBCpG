.. CluBCpG documentation master file, created by
   sphinx-quickstart on Fri May 17 09:19:49 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================================================
CluBCpG: Cluster-Based analysis of CpG methylation
===================================================


What is CluBCpG?
===================================
CluBCpG is a software package built to analyze whole genome bisulfite sequencing (WGBS) data.
This toolkit will divide each chromosome into small user-defined intervals, extract all WGBS reads within those intervals,
cluster them based on identity, and write a final report to the use containing all identified CpG methylation patterns.

CluBCpG supports both paired-end and single-end data and any read-length. However, a read length of at least
100 base-pairs is recommended

CluBCpG is written 100% in Python.

CluBCpG was created to be primarily utilized as command line based tools, but it does provide a few APIs which may be
useful on their own for researchers. Those APIs are documented in :doc:`API`.


.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: User Guide Contents:

   intro.rst
   usage.rst

.. toctree::
   :maxdepth: 2
   :caption: Developer Reference:

   API.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
