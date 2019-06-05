========================================
Understanding the CluBCpG output data
========================================

Coverage output
================
The output from ``clubcpg-coverage`` looks like this:

.. list-table:: clubcpg-coverage output
    :widths: 10 5 5
    :header-rows: 1

    * - bin_id
      - n_reads
      - n_cpgs
    * - chr19_61300
      - 4
      - 1
    * - chr19_61400
      - 4
      - 2
    * - chr19_89800
      - 25
      - 4

Each row represents 1 **bin**.

.. NOTE::
    A header row is shown on this table for clarity, however the real csv file generated will not have a header row.

* bin_id
    Represents the unique bin in the genome. The underscore character (``__``) separates the chromosome and the genomic
    coordinate. The genomic coordiate represents the end-point of a bin.

        * ex: If ``--bin_size`` was set to 100, chr19_61300 would represent chr19:61200-61300

* n_reads
    Number of reads which fully cover **all** CpGs within the bin

* n_cpgs
    Number of CpGs within the bin


Cluster output
================
The output from ``clubcpg-cluster`` looks like this:

.. list-table:: clubcpg-cluster output
    :widths: 10 10 10 10 10 10 10 10
    :header-rows: 1

    * - bin_id
      - input_label
      - methylation
      - class_label
      - read_number
      - cpg_number
      - cpg_pattern
      - class_split
    * - chr2_10700
      - AB
      - 0.8333
      - 0
      - 6
      - 6
      - 1;1;1;1;1;0
      - A=5,B=1
    * - chr2_10700
      - A
      - 1
      - 1
      - 12
      - 6
      - 1;1;1;1;1;1
      - A=12

Each row represents 1 **cluster**.

* bin_id
    Represents the unique bin in the genome. The underscore character (``__``) separates the chromosome and the genomic
    coordinate. The genomic coordiate represents the end-point of a bin.

* input_label
    Represents the input file this cluster was found in

    * Single-file mode
        * This will show the BAM file name specified with the ``-a`` flag.

    * Two-file mode:
        * A = Input BAM specified in the ``-a`` flag
        * B = Input BAM specified in the ``-b``
        * AB = found in both input BAM files

* methylation
    The methylation level of the CpG pattern found in this cluster

* class_label
    A unique identifier for each cluster found within a given bin (only unique within a bin)

* read_number
    The number of reads found within this cluster

* cpg_number
    The number of CpGs in this cluster

* cpg_pattern
    The methylation pattern of the cluster

        * 1 = methylated
        * 0 = unmethylated

* class_split
    The breakdown of how many reads in a given cluster came from each input BAM.

        * Useful for AB clusters.
        * If cluster is A or B, this number will match the read_number column.

.. WARNING::
    If clustering was performed with ``--remove_noise False`` (this is ``True`` by deafult) you may find clusters with a
    class_label of -1.

    These represents noise reads, which are reads where a given CpG patterns was only observed **ONCE**.

    The *read_number* column may be of interest to you, but the other columns do **NOT** represent true values. They will
    only represent one of the noise patterns found. These values may be set to ``null`` in future versions of CluBCpG.

