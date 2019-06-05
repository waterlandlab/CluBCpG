===========================
Using CluBCpG with PReLIM
===========================

Introduction
=============
PReLIM exists as its own stand-alone package, but for simplicity and compatibility, a version of PReLIM comes
bundled with CluBCpG.

Imputation with PReLIM is taken care of behind the scenes and CluBCpG includes three command line scripts to perform
analysis with imputation.

Usage is almost identical to standard CluBCpG usage, but includes one extra step and requires the use of a
couple extra command line flags.

The command-tools provided are ``clubcpg-impute-coverage``, ``clubcpg-impute-train``, and ``clubcpg-impute-cluster``

Typical imputation workflow
============================

Calculate bin coverage
***********************

Use ``clubcpg-coverage`` to calculate bin coverage as performed in the standard workflow. (:ref:`typical_workflow_label`)

Filter outputs
***************

Filter these csv outputs for >= 1 reads and >= 2 CpGs. The following one-liner will filter it correctly:

.. code-block:: bash

        cat CompleteBins.yourfilename.chr19.csv | awk -F "," '$2>=1 && $3>=2' > CompleteBins.yourfilename.chr19.filtered.csv

But as before, you can filter this with any other method you like.

.. NOTE::
    PReLIM requires at least 1 read fully covering all CpGs in a bin


