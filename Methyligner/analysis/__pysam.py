##/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## Python shit
import os
import pysam
import shutil
import itertools
import pysamstats
import subprocess
import numpy as np
import logging as log
import multiprocessing

## Backend shit
THREADS = multiprocessing.cpu_count()
from ..__backend import force_mkdir
from ..__backend import Colour as clr

class Quantification:
    def __init__(self, sequencepair_object, instance_params):
    	self.sequencepair_object = sequencepair_object
    	self.instance_params = instance_params
    	self.reference_sequence = self.sequencepair_object.get_referenceidx()
    	self.forward_aln = self.sequencepair_object.get_forward_convassembly()
    	self.reverse_aln = self.sequencepair_object.get_reverse_convassembly()
    	self.sample_root = sequencepair_object.get_label()
    	self.target_output = sequencepair_object.get_analysispath()

    	## run functions
    	for assembly in [(self.forward_aln, 'R1'), (self.reverse_aln, 'R2')]:
            self.determine_request(assembly[0], assembly[1])

    def determine_request(self, bamfi, orientation):

    	"""
    	Determine which aspects of pysamstats the user wishes to be run
    	and on how many reads to execute said stages
    	"""
    	quantify_variation = self.instance_params.config_dict['analysis_flags']['@quant_variation']
    	quantify_mapq = self.instance_params.config_dict['analysis_flags']['@quant_mapq']
    	quantify_baseq = self.instance_params.config_dict['analysis_flags']['@quant_baseq']
    	read_depth_limit = int(self.instance_params.config_dict['analysis_flags']['@read_depth_limit'])

    	## get read count as depth limit if user specified 0
    	## else use user specified depth limit
    	utilised_limit = 0
    	if read_depth_limit == 0:
    		utilised_limit = int(pysam.idxstats(bamfi).split('\n')[0].split('\t')[2])
    	else:
    		utilised_limit = read_depth_limit

    	# determine methylation
    	methylation_data = self.pysamstats_scraper(bamfi, orientation, '_methylationReport.txt', 'coverage_gc', 3, utilised_limit)
    	## assign objects for this orientation
    	if orientation == 'R1': self.sequencepair_object.set_forward_methylation(methylation_data)
    	if orientation == 'R2': self.sequencepair_object.set_reverse_methylation(methylation_data)

    	## determine variation
    	if quantify_variation: 
    		variation_data = self.pysamstats_scraper(bamfi, orientation, '_variationReport.txt', 'variation', 21, utilised_limit)
    		## assign objects for this orientation
    		if orientation == 'R1': self.sequencepair_object.set_forward_variation(variation_data)
    		if orientation == 'R2': self.sequencepair_object.set_reverse_variation(variation_data)

    	## determine mapQ
    	if quantify_baseq:
    		baseq_data = self.pysamstats_scraper(bamfi, orientation, '_baseqReport.txt', 'baseq', 4, utilised_limit)
    		## assign objects for this orientation
    		if orientation == 'R1': self.sequencepair_object.set_forward_baseq(baseq_data)
    		if orientation == 'R2': self.sequencepair_object.set_reverse_baseq(baseq_data)

   		## determine baseQ
   		if quantify_mapq:
			mapq_data = self.pysamstats_scraper(bamfi, orientation, '_mapqReport.txt', 'mapq', 8, utilised_limit)
			## assign objects for this orientation
			if orientation == 'R1': self.sequencepair_object.set_forward_mapq(mapq_data)
			if orientation == 'R2': self.sequencepair_object.set_reverse_mapq(mapq_data)

    def pysamstats_scraper(self, bamfi, orientation, outfi_string, requested_analysis, blank_range, read_depth):

    	##
    	## Generate data for the current desired function
    	## requested_analysis == specific stats the user wants to run
    	utilised_reference = self.sequencepair_object.get_referencefile()
    	target_output_file = os.path.join(self.target_output, '{}{}'.format(orientation, outfi_string))
    	target_output_object = open(target_output_file, 'w')
    	target_subprocess = subprocess.Popen(['pysamstats', '--type', requested_analysis, '--max-depth',
    	 str(read_depth), '--fasta', utilised_reference, bamfi],
    	 stdout=target_output_object, stderr=subprocess.PIPE)
    	subprocess_stderr = target_subprocess.communicate()[1]; target_subprocess.wait(); target_output_object.close()

    	##
    	## Filter results for the CPG region's relevant positions
    	target_region = self.instance_params.config_dict['MethRegion'][0]
    	target_positions = self.sequencepair_object.methylation_regions(target_region)
    	np.warnings.filterwarnings('ignore') ##ignore numpycomplaining about empty file if pysamstats returned null
    	analysis_data = np.ndarray.tolist(np.genfromtxt(target_output_file, delimiter='\t', dtype=None, encoding='utf8', skip_header=1))
    	if len(analysis_data) == 0:
    		##empty file, no analysis reported
    		for position in target_positions:
    			header = ['{}{}'.format(target_region, 'NullAnalysisReported'), position]; footer = ['0']*blank_range
    			analysis_data.append(header+footer)
    	else:
    		## file had analysis
    		## remove unwanted positions for this CPG region, iterate backwards
    		for i in xrange(len(analysis_data)-1, -1, -1):
    			element = analysis_data[i]
    			if element[1] not in target_positions:
    				del analysis_data[i]

    		## check leftover analysis positions for missing values
    		for position, result in itertools.izip_longest(target_positions, analysis_data):
    			try:
    				if position != result[1]:
    					mismatch_idx = analysis_data.index(result)
    					header = ['{}{}'.format(target_region, 'NullAnalysisReported'), position]; footer = ['0']*blank_range
    					analysis_data.insert(mismatch_idx, header+footer)
    			except TypeError:
    				## variation data had a None within iteration of valid positions
    				header = ['{}{}'.format(target_region, 'NullAnalysisReported'), position]; footer = ['0']*blank_range
    				analysis_data.append(header+footer)

    	return analysis_data