from __future__ import print_function

"""
Author: 
Jack Duryea
Waterland Lab
Computational Epigenetics Section
Baylor College of Medicine

April 2018

A Python object representing information about a 
bin and its meta data. Very useful for storing information and working with it.
Can be saved as a pickle file.
"""
import numpy as np

class Bin():
	""" 
	Constructor for a bin

	Inputs:

		matrix: numpy array, the bin's CpG matrix.
		binStartInc: integer, the starting, inclusive, chromosomal index of the bin.
		binEndInc: integer, the ending, inclusive, chromosomal index of the bin.
		cpgPositions: array of integers, the chromosomal positions of the CpGs in the bin.
		sequence: string, nucleotide sequence (A,C,G,T)
		encoding: array, a reduced representation of the bin's CpG matrix
		missingToken: integer, the token that represents missing data in the matrix.
		chromosome: string, the chromosome this bin resides in.
		binSize: integer, the number of base pairs this bin covers
		species: string, the speices this bin belongs too.
		verbose: boolean, print warnings, set to "false" for no error checking and faster speed

		tag1: anything, for custom use.
		tag2: anything, for custom use.
	"""
	def __init__(self, 
			matrix, 
			binStartInc, 
			binEndInc, 
			cpgPositions, 
			sequence="",
			encoding=None, 
			missingToken= -1, 
			chromosome="19", 
			binSize=100, 
			species="HG38", 
			verbose=True, 
			tag1=None, 
			tag2=None):


		self.cpgDensity = matrix.shape[1]
		self.readDepth = matrix.shape[0]
		
		if verbose:
			assert binSize > 0, "invalid bin size"
			assert binStartInc < binEndInc, "invalid start and end indices"
			assert binEndInc - binStartInc == binSize - 1
			# TODO: add more assertions
			assert len(cpgPositions) == self.cpgDensity, "wrong number of positions"

			if (not (species == "HG38")) and (not (species == "MM10")):
				print("Warning, you are not supplying a common species type. You've been warned")

			assert self.readDepth > 0, "invalid read depth"

		self.matrix = matrix
		self.binStartInc = binStartInc
		self.binEndInc = binEndInc
		self.cpgPositions = cpgPositions
		self.sequence = sequence
		self.missingToken = missingToken
		self.chromosome = chromosome
		self.binSize = binSize
		self.species = species
		self.tag1 = tag1
		self.tag2 = tag2


	






