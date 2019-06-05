===========================
Using CluBCpG with PReLIM
===========================

Introduction
=============
PReLIM exists as its own stand-alone package, but for simplicity and compatibility, a version of PReLIM comes
bundled with CluBCpG.

Imputation with PReLIM is taken care of behind the scenes and CluBCpG includes three command line scripts to perform
analysis with imputation.

Usage is almost identical to standard CluBCpG pipeline, but includes one extra step and requires the use of a
couple extra command line flags.

The documentation included in this section mostly highlights the differences and additional steps needed to run CluBCpG with
PReLIM imputation.

.. NOTE::
    It is highly recommended that you have read :ref:`using_CluBCpG_label` first.

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

.. NOTE::
    Differing from the typical workflow, if you intend to analyze two BAM files with imputation, you will want to perform this process on both BAM files.
    This is because you will train an imputation model using both BAM files separately. This improves the accuracy of the imputations.

Train a PReLIM model
*********************

Now ``clubcpg-impute-train`` can be used to automatically train multiple imputation models using PReLIM. For each BAM file, a random set of bins
will be selected from the filtered coverage file provided.

The number of randomly sampled bins can be set using the ``-l`` flag.

.. WARNING::
    Increasing the number of bins will drastically increase the run-time of this process and can consume a lot of memory, potentially crashing PReLIM.
    When testing we found there was no improvement in accuracy above 10,000 bins.

Models bins containing 2-5 CpGs will be automatically trained and saved in the folder designated by ``-o``. Do **NOT** rename them.

Compute coverage gains
************************

Now with the trained models saved, the post-imputation coverage can be calculated using ``clubcpg-impute-coverage``. This
functions almost identically to ``clubcpg-coverage`` except:

* you need to provide the folder containing the saved models using the ``-m`` flag.

* You also need to provide the coverage file filtered for >= 1 reads with the ``-c`` flag.

    * This accelerates the process by skipping bins which cannot be imputed upon do to lack of coverage.

This will create an output file identical to ``clubcpg-coverage``, except the values represent the number of reads post-imputation

Filter imputed coverage output
*******************************

Filter this csv output however you wish. Again, we recommend 10 reads and 2 cpgs. See :ref:`typical_filter_label` for more details.

Perform clustering with imputation
***********************************

Using ``clubcpg-impute-cluster`` you can now perform read clustering. This also functions almost identically to the standard method.
However, you must also point it to the folder containing the models.

Here you use the ``--models_A`` and ``--models_B`` flags to point the tool to PReLIM models for input ``-a`` and ``-b`` respectively. Imputation on each file
is performed independent of the other.

Command line tools
===================
These options can also be viewed by running ``--help`` after each tool on the command line.

.. autoprogram:: clubcpg-impute-train:arg_parser
    :prog: clubcpg-impute-train


.. autoprogram:: clubcpg-impute-coverage:arg_parser
    :prog: clugcpg-impute-coverage


.. autoprogram:: clubcpg-impute-cluster:arg_parser
    :prog: clubcpg-impute-cluster



