import pysam
import argparse
import copy
import pandas as pd
import numpy as np
import sys
from collections import defaultdict


class BamFileReadParser():

    def __init__(self, bamfile, quality_score, read1_5=None, read1_3=None,
                 read2_5=None, read2_3=None, no_overlap=True):

        self.mapping_quality = quality_score
        self.bamfile = bamfile
        self.read1_5 = read1_5
        self.read1_3 = read1_3
        self.read2_5 = read2_5
        self.read2_3 = read2_3
        self.full_reads = []
        self.read_cpgs = []
        self.no_overlap = no_overlap

        if read1_5 or read2_5 or read1_3 or read2_3:
            self.mbias_filtering = True
        else:
            self.mbias_filtering = False

        self.OpenBamFile = pysam.AlignmentFile(bamfile, 'rb')
        # Check for presnece of index file
        assert self.OpenBamFile.check_index(), "Can't find index file. Please run samtools index to generate it."

    # Get reads from the bam file, extract methylation state
    def parse_reads(self, chromosome, start, stop):
        reads = []
        for read in self.OpenBamFile.fetch(chromosome, start, stop):
            if read.mapping_quality >= self.mapping_quality:
                reads.append(read)

        read_cpgs = []

        for read in reads:
            reduced_read = []
            # Join EVERY XM tag with its aligned_pair location
            for pair, tag in zip(read.get_aligned_pairs(), read.get_tag('XM')):
                if pair[1]:
                    if read.flag == 83 or read.flag == 163:
                        reduced_read.append((pair[1] - 1, tag))
                    else:
                        reduced_read.append((pair[1], tag))
                else:
                    continue

            # if MBIAS was set, slice the joined list
            if self.mbias_filtering:
                if read.is_read1:
                    mbias_5_prime = self.read1_5
                    # note taking the NEGATIVE of the value for the 3-prime
                    mbias_3_prime = -self.read1_3
                    if mbias_3_prime == 0:
                        mbias_3_prime = None
                    reduced_read = reduced_read[mbias_5_prime:mbias_3_prime]
                if read.is_read2:
                    mbias_5_prime = self.read2_5
                    mbias_3_prime = -self.read2_3
                    if mbias_3_prime == 0:
                        mbias_3_prime = None
                    reduced_read = reduced_read[mbias_5_prime:mbias_3_prime]

            read_cpgs.append(reduced_read)

        self.full_reads = reads
        self.read_cpgs = read_cpgs

        # Correct overlapping paired reads if set, this is default behavior
        if self.no_overlap:
            try:
                read_cpgs = self.fix_read_overlap(reads, read_cpgs)
            except AttributeError:
                print("Could not determine read 1 or 2. {}:{}-{}".format(chromosome, start, stop))
                sys.stdout.flush()

        # Filter the list for positions between start-stop and CpG (Z/z) tags
        output = []
        for read_cpg in read_cpgs:
            temp = []
            for pos, tag in read_cpg:
                if pos and pos >= start and pos <= stop and (tag == 'Z' or tag == 'z'):
                    temp.append((pos, tag))
            output.append(temp)

        return output

    def create_matrix(self, read_cpgs):
        series = []
        for data in read_cpgs:
            positions = []
            statues = []
            for pos, status in data:
                positions.append(pos)
                statues.append(status)
            series.append(pd.Series(statues, positions))
        matrix = pd.concat(series, axis=1)
        matrix = matrix.replace('Z', 1)
        matrix = matrix.replace('z', 0)

        return matrix.T

    def fix_read_overlap(self, full_reads, read_cpgs):
        """
        Takes pysam reads and read_cpgs generated during parse reads and removes any
        overlap between read1 and read2. If possible it also stitches read1 and read2 together to create
        a super read.
        :param full_reads:
        :param read_cpgs:
        :return: A list in the same format as read_cpgs input, but corrected for paired read overlap
        """
        # data for return
        fixed_read_cpgs = []
        # Combine raw reads and extracted tags
        combined = []
        for read, state in zip(full_reads, read_cpgs):
            combined.append((read, state))

        # Get names of all the reads present
        query_names = [x.query_name for x in full_reads]

        # Match paired reads by query_name
        tally = defaultdict(list)
        for i, item in enumerate(query_names):
            tally[item].append(i)

        for key, value in sorted(tally.items()):
            # A pair exists, process it
            if len(value) == 2:
                # Set read1 and read2 correctly
                if combined[value[0]][0].is_read1:
                    read1 = combined[value[0]]
                    read2 = combined[value[1]]

                elif combined[value[1]][0].is_read1:
                    read1 = combined[value[1]]
                    read2 = combined[value[0]]

                # both reads have same value, this shouldn't be. Drop one completely, dont
                # bother with overlap
                elif combined[value[0]][0].is_read1 == combined[value[1]][0].is_read1:
                    fixed_read_cpgs.append(combined[value[0]][1])
                    continue


                else:
                    raise AttributeError("Could not determine read 1 or read 2")

                # Find amount of overlap
                amount_overlap = 0
                r1_bps = [x[0] for x in read1[1]]
                r2_bps = [x[0] for x in read2[1]]

                if min(r1_bps) < min(r2_bps):
                    trim_direction = 5
                    for bp in r2_bps:
                        if bp and bp in r1_bps:
                            amount_overlap += 1
                else:
                    trim_direction = 3
                    for bp in r1_bps:
                        if bp and bp in r2_bps:
                            amount_overlap += 1


                # remove the overlap by trimming or discarding
                if amount_overlap == len(read2[1]):
                    # discard read 2, only append read 1
                    fixed_read_cpgs.append(read1[1])
                else:
                    # trim overlap
                    if trim_direction == 5:
                        new_read2_cpgs = read2[1][amount_overlap:]
                    elif trim_direction == 3:
                        new_read2_cpgs = read2[1][:-amount_overlap]
                    # stitch together read1 and read2
                    read1[1].extend(new_read2_cpgs)
                    fixed_read_cpgs.append(read1[1])

            elif len(value) == 1:
                # No pair, add to output
                fixed_read_cpgs.append(combined[value[0]][1])

        return fixed_read_cpgs


