from __future__ import print_function

"""
Author: 
Jack Duryea
Waterland Lab
Computational Epigenetics Section
Baylor College of Medicine

Created April 2018

Updated April 11 2019: use random forests as model


PReLIM: Preceise Read Level Imputation of Methylation

PReLIM imputes missing CpG methylation
states in CpG matrices.

"""

# standard imports
import numpy as np
from tqdm import tqdm
import copy
import time

from collections import defaultdict
import random

# sklearn imports
from sklearn.preprocessing import normalize
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV

# Pickle
try:
	import cPickle as p
except ImportError:
	import pickle as p


# TODO: most of these fields are redundant in our application
class CpGBin():
	""" 
	Constructor for a bin

	"""
	def __init__(self, 
			matrix, 
			#relative_positions
			binStartInc=None, 
			binEndInc=None, 
			cpgPositions=None, 
			sequence="",
			encoding=None, 
			missingToken= -1, 
			chromosome=None, 
			binSize=100, 
			species="MM10", 
			verbose=True, 
			tag1=None, 
			tag2=None):
		"""

		:param matrix: numpy array, the bin's CpG matrix.
		:param binStartInc: integer, the starting, inclusive, chromosomal index of the bin.
		:param binEndInc: integer, the ending, inclusive, chromosomal index of the bin.
		:param cpgPositions: array of integers, the chromosomal positions of the CpGs in the bin.
		:param sequence: string, nucleotide sequence (A,C,G,T)
		:param encoding: array, a reduced representation of the bin's CpG matrix
		:param missingToken: integer, the token that represents missing data in the matrix.
		:param chromosome: string, the chromosome this bin resides in.
		:param binSize: integer, the number of base pairs this bin covers
		:param species: string, the speices this bin belongs too.
		:param verbose: boolean, print warnings, set to "false" for no error checking and faster speed
		:param tag1: anything, for custom use.
		:param tag2: anything, for custom use.

		"""


		self.cpgDensity = matrix.shape[1]
		self.readDepth = matrix.shape[0]
		
		

		self.matrix = np.array(matrix, dtype=float)
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


class PReLIM():
	"""
	PReLIM imputation class to handle training and predicting from models.

	"""
	def __init__(self, cpgDensity=2):
		self.model = None
		self.cpgDensity = cpgDensity
		self.METHYLATED = 1
		self.UNMETHYLATED = 0
		self.MISSING = -1
		self.methylated = 1
		self.unmethylated = 0
		self.unknown = -1

	# Train a model
	def train(self, bin_matrices, model_file="no", verbose=False):
		"""
		bin_matrices: list of cpg matrices

		model_file, string,      The name of the file to save the model to. 
			If None, then create a file name that includes a timestamp.
			If you don't want to save a file, set this to "no"
		"""
		# bin_matrices: a list of cpg matrices 
		X,y = self.get_X_y(bin_matrices, model_file=model_file, verbose=False)
		
		# Train the neural network model
		self.fit(X,y, model_file=model_file, verbose=verbose)

	def fit(self,
			X_train,
			y_train,
			n_estimators = [10, 50, 100, 500, 1000],
			cores = -1,
			max_depths = [1, 5, 10, 20, 30],
			model_file=None,
			verbose=False
			):
		"""
		Inputs: 
		1. X_train,     numpy array, Contains feature vectors.
		2. y_train,     numpy array, Contains labels for training data.
		3. n_estimators, list, the number of estimators to try during a grid search.
		4. max_depths, list, the maximum depths of trees to try during a grid search.
		5. cores, the number of cores to use during training, helpful for grid search.
		6. model_file, string,      The name of the file to save the model to. 
			If None, then create a file name that includes a timestamp.
			If you don't want to save a file, set this to "no"
	
		5-fold validation is built into the grid search

		Outputs: 
		The trained model

		Usage: 
		model.fit(X_train, y_train)	
		"""


		grid_param = {  
		 "n_estimators": n_estimators,
		  "max_depth": max_depths,
		}

		# Note: let the grid search use a lot of cores, but only use 1 for each forest
		# since dispatching can take a lot of time
		rf = RandomForestClassifier(n_jobs=1)
		self.model = GridSearchCV(rf, grid_param, n_jobs=cores, cv=5, verbose=verbose)
		self.model.fit(X_train, y_train)


		# save the model

		if model_file == "no":
			return self.model

		if not model_file:
			model_file = "PReLIM_model" + str(time.time())

		p.dump(self.model, open(model_file,"wb"))

		return self.model




	# Feature collection directly from bins
	def get_X_y(self, bin_matrices, model_file=None, verbose=False):
		bins = []

		# convert to bin objects for ease of use
		for matrix in bin_matrices:
			mybin = CpGBin( matrix=matrix )
			bins.append( mybin )
		
		# find bins with no missing data
		complete_bins = _filter_missing_data( bins )
		random.shuffle( complete_bins )
		
		# apply masks
		masked_bins = _apply_masks( complete_bins, bins )

		# extract features
		X, y = self._collectFeatures( masked_bins ) 
		return X, y 


	# Return a vector of predicted classes 
	def predict_classes(self, X):
		"""
		Inputs: 
		1. X, numpy array, contains feature vectors
		
		Outputs: 
		1. 1-d numpy array of prediction values

		Usage: 
		y_pred = CpGNet.predict_classes(X)  

		"""
		return self.model.predict(X)
	
	# Return a vector of probabilities for methylation
	def predict(self, X):
		"""
		Inputs: 
		1. X, numpy array, contains feature vectors
		
		Outputs: 
		1. 1-d numpy array of predicted class labels

		Usage: 
		y_pred = CpGNet.predict(X)  

		"""
		return self.model.predict_proba(X)[:,1]


	def predict_proba(self, X):
		"""
		Inputs: 
		1. X, numpy array, contains feature vectors
		
		Outputs: 
		1. 1-d numpy array of class predictions 

		Usage: 
		y_pred = CpGNet.predict(X)  

		"""
		return self.model.predict_proba(X)[:1]


	# Load a saved model
	def loadWeights(self, model_file):
		"""
		Inputs:
		1. model_file, string, name of file with a saved model

		Outputs:
		None

		Effects:
		self.model is loaded with the provided weights
		"""
		self.model = p.load(open(model_file,"rb"))

	




	def _get_imputation_features(self,matrix):
		'''
		Returns a vector of features needed for the imputation of this matrix

		Inputs: 
		1. matrix, a 2d np array, dtype=float, representing a CpG matrix, 1=methylated, 0=unmethylated, -1=unknown

		Outputs:
		1. A feature vector for the matrix
		'''

		X = []

		numReads = matrix.shape[0]
		density = matrix.shape[1]

		nan_copy = np.copy(matrix)
		nan_copy[nan_copy == -1] = np.nan
		column_means = np.nanmean(nan_copy, axis=0)
		row_means = np.nanmean(nan_copy, axis=1)
		
		encoding = self._encode_input_matrix(matrix)[0]

		for i in range(numReads):
			for j in range(density):
				observed_state = matrix[i, j]
				
				if observed_state != -1:
					continue

				row_mean = row_means[i]
				col_mean = column_means[j]
				
				row = np.copy(matrix[i])
				row[j] = -1

				data = [row_mean] + [col_mean] +  [i, j] + list(row) +  list(encoding)
				X.append(data)

		X = np.array(X)

		return X

	# Imputes missing values in Bins
	def impute(self, matrix):
		"""
		Inputs: 
		1. matrix, a 2d np array, dtype=float, representing a CpG matrix, 1=methylated, 0=unmethylated, -1=unknown
		
		Outputs: 
		1. A 2d numpy array with predicted probabilities of methylation

		"""

		X = self._get_imputation_features(matrix)
		
		if len(X) == 0: # nothing to impute
			return matrix

		predictions = self.predict(X)

		k = 0 # keep track of prediction index for missing states
		predicted_matrix = np.copy(matrix)
		for i in range(predicted_matrix.shape[0]):
			for j in range(predicted_matrix.shape[1]):
				if predicted_matrix[i, j] == -1:
					predicted_matrix[i, j] = predictions[k]
					k += 1

		return predicted_matrix







	def impute_many(self, matrices):
		'''
		Imputes a bunch of matrices at the same time to help speed up imputation time.

		Inputs:

		1. matrices: array-like (i.e. list), where each element is
		a 2d np array, dtype=float, representing a CpG matrix, 1=methylated, 0=unmethylated, -1=unknown

		Outputs:

		1. A List of 2d numpy arrays with predicted probabilities of methylation for unknown values.
		'''

		# Extract all features for all matrices so we can predict in bulk, this is where the speedup comes from
		
		X = np.array([features for matrix_features in [self._get_imputation_features(matrix) for matrix in matrices] for features in matrix_features])
		
		if len(X) == 0:
			return matrices
		
		predictions = self.predict(X)


		predicted_matrices = []

		# TODO: lots of for-loops here, could be sped up?

		k = 0 # keep track of prediction index for missing states, order is crucial!
		for matrix in matrices:
			predicted_matrix = np.copy(matrix)
			for i in range(predicted_matrix.shape[0]):
				for j in range(predicted_matrix.shape[1]):
					if predicted_matrix[i, j] == -1:
						predicted_matrix[i, j] = predictions[k]
						k += 1
			predicted_matrices.append(predicted_matrix)

		return predicted_matrices







	### Helper functions, for private use only ###

	# Returns a matrix encoding of a CpG matrix
	def _encode_input_matrix(self, m):
		matrix = np.copy(m)
		n_cpgs = matrix.shape[1]
		matrix += 1  # deal with -1s
		base_3_vec = np.power(3, np.arange(n_cpgs - 1, -1, -1))
		#
		encodings = np.dot(base_3_vec, matrix.T)
		#
		encoded_vector_dim = np.power(3, n_cpgs)
		encoded_vector = np.zeros(encoded_vector_dim)
		#
		for x in encodings:
			encoded_vector[int(x)] += 1
		#
		num_reads = encodings.shape[0]
		#
		# Now we normalize
		encoded_vector_norm = normalize([encoded_vector], norm="l1")
		return encoded_vector_norm[0], num_reads

	# finds the majority class of the given column, discounting the current cpg
	
	def _get_column_mean(self, matrix, col_i, current_cpg_state):
		sub = matrix[:, col_i]
		return self._get_mean(sub, current_cpg_state)

	# finds the majority class of the given read, discounting the current cpg
	def _get_read_mean(self, matrix, read_i, current_cpg_state):
		sub = matrix[read_i, :]
		return self._get_mean(sub, current_cpg_state)

	def _get_mean(self, sub_matrix, current_cpg_state):
		num_methy = np.count_nonzero(sub_matrix == self.METHYLATED)
		num_unmethy = np.count_nonzero(sub_matrix == self.UNMETHYLATED)

		if current_cpg_state == self.METHYLATED:
			num_methy -= 1
		num_methy = max(0, num_methy)
		if current_cpg_state == self.UNMETHYLATED:
			num_unmethy -= 1
		num_unmethy = max(0, num_unmethy)
		if float(num_methy + num_unmethy) == 0:
			return -2

		return float(num_methy) / float(num_methy + num_unmethy)

	# Returns X, y
	# note: y can contain the labels 1,0, -1
	def _collectFeatures(self, bins):
		X = []
		Y = []
		for Bin in tqdm(bins):
			observed_matrix = Bin.tag2["observed"]
			truth_matrix = Bin.tag2["truth"]
			encoding = self._encode_input_matrix(observed_matrix)[0]

			numReads = observed_matrix.shape[0]
			density = observed_matrix.shape[1]
			#positions = Bin.cpgPositions
			nan_copy = np.copy(observed_matrix)
			

			nan_copy[nan_copy == -1] = np.nan 
			column_means = np.nanmean(nan_copy,axis=0)
			row_means = np.nanmean(nan_copy,axis=1)

			for i in range(numReads):
				for j in range(density):
					observed_state = observed_matrix[i,j]
					if observed_state != -1:
						continue

					state = truth_matrix[i,j]
					Y.append(state)
					
					
					
					row_mean = row_means[i]
					col_mean = column_means[j]
					# j is the current index in the row
					
					# encoding is the matrix encoding vector
					# differences is the difference in positions of the cpgs
					row = np.copy(observed_matrix[i])
					row[j] = -1

					data = [row_mean] + [col_mean] +  [i, j] + list(row) +  list(encoding)
					X.append(data)


		X = np.array(X)
		Y = np.array(Y)
		Y.astype(int)
		return X, Y



# Helper functions

# returns a list of bins similar to the input
# but matrix rows with missing values are removed
def _filter_bad_reads(bins):
	filtered_bins = []
	for Bin in bins:
		newBin = copy.deepcopy(Bin)
		matrix = newBin.matrix

		# find rows with missing values
		counts = np.count_nonzero(matrix == -1, axis=1)
		idx = counts == 0
		matrix_filtered = matrix[idx]
		newBin.matrix = matrix_filtered
		filtered_bins.append(newBin)

	return filtered_bins

# returns a mapping of dimensions to list of masks that can be used on data
# of that size.
# the missing pattern is in matrix form.
# -1 is missing, 2 is known
def _extract_masks( bins):
	masks = defaultdict(lambda: [])
	for Bin in tqdm(bins):
		matrix = np.copy(Bin.matrix)
		matrix[matrix >= 0] = 2

		#min_missing = 10
		min_missing = 1 # must have at least 1 missing value
		if np.count_nonzero(matrix == -1) >= min_missing:
			masks[matrix.shape].append(matrix)

	return masks

  
def _apply_masks( filtered_bins, all_bins ):
	
	masks = _extract_masks( all_bins )
	ready_bins = []

	for Bin in filtered_bins:
		truth_matrix = Bin.matrix
		m_shape = truth_matrix.shape
		if m_shape in masks:
			if len( masks [ m_shape ] ) > 0:
				mask = random.choice(masks[m_shape])
				observed = np.minimum(truth_matrix, mask)
				Bin.tag2 = {"truth":truth_matrix, "observed":observed, "mask":mask}
				ready_bins.append(Bin)

	return ready_bins

# get a set of bins with no missing data
def _filter_missing_data( bins, min_read_depth=1 ):
	cpg_bins_complete = _filter_bad_reads(bins)
	# secondary depth filter
	cpg_bins_complete_depth = [bin_ for bin_ in cpg_bins_complete if bin_.matrix.shape[0] >= min_read_depth]
	return cpg_bins_complete_depth


