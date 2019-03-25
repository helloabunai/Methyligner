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

        self.determine_methylation(bamfi, orientation)
    	if quantify_variation: self.determine_variation(bamfi, orientation)
    	if quantify_mapq: self.determine_mapq(bamfi, orientation)
    	if quantify_baseq: self.determine_baseq(bamfi, orientation)

    def determine_variation(self, bamfi, orientation):

        ##
        ## Generate variation table
        utilised_reference = self.sequencepair_object.get_referencefile()
        variation_output = os.path.join(self.target_output, '{}_variation_report.txt'.format(orientation))
        variation_file = open(variation_output, 'w')
        variation_process = subprocess.Popen(['pysamstats', '--type', 'variation', '--fasta', utilised_reference, bamfi], stdout=variation_file)
        variation_stderr = variation_process.communicate()[1]; variation_process.wait(); variation_file.close()

        ##
        ## Filter for current region's position(s)
        target_region = self.instance_params.config_dict['MethRegion'][0]
        methregion_positions = self.sequencepair_object.methylation_regions(target_region)
        np.warnings.filterwarnings('ignore') ## ignore numpy complaining about empty file if we skip header on variation_report with no results
        variation_data = np.ndarray.tolist(np.genfromtxt(variation_output, delimiter='\t', dtype=None, encoding="utf8", skip_header=1))
        ## HEADER:: chrom, pos, ref, reads_all, reads_pp, matches, matches_pp, mismatches, mismatches_pp, deletions, deletions_pp, insertions, insertions_pp, A, A_pp, C, C_pp, T, T_pp, G, G_pp, N, N_pp
        ## Filter all processed positions to those relevant to current methylation region
        if len(variation_data) == 0:
            ## empty file, no variation reported for this alignment
            for position in methregion_positions:
                temp1 = ['{}{}'.format(target_region, 'NullVariationReported'), position]; temp2 = ['0']*21
                variation_data.append(temp1+temp2)
        else:
            ## file had variation report
            ## remove positions from source data not in the current CPG region
            ## iterate BACKWARDS to remove faster
            for i in xrange(len(variation_data) - 1, -1, -1):
                element = variation_data[i]
                if element[1] not in methregion_positions:
                    del variation_data[i]

            ## check for positions missing for current CPG region data
            for position, result in itertools.izip_longest(methregion_positions, variation_data):
                try:
                    if position != result[1]:
                        mismatch_idx = variation_data.index(result)
                        temp1 = ['{}{}'.format(target_region, 'NullVariationReported'), position]; temp2 = ['0']*21
                        variation_data.insert(mismatch_idx, temp1+temp2)
                except TypeError:
                    ## variation_data had a None within iteration
                    temp1 = ['{}{}'.format(target_region, 'NullVariationReported'), position]; temp2 = ['0']*21
                    variation_data.append(temp1+temp2)

        ## assign to object
        if orientation == 'R1': self.sequencepair_object.set_forward_variation(variation_data)
        if orientation == 'R2': self.sequencepair_object.set_reverse_variation(variation_data)

    def determine_methylation(self, bamfi, orientation):

        ##
        ## Generate methylation data
        utilised_reference = self.sequencepair_object.get_referencefile()
        methylation_output = os.path.join(self.target_output, '{}_methylation_report.txt'.format(orientation))
        methylation_file = open(methylation_output, 'w')
        methylation_process = subprocess.Popen(['pysamstats', '--type', 'coverage_gc', '--fasta', utilised_reference, bamfi], stdout=methylation_file)
        methylation_stderr = methylation_process.communicate()[1]; methylation_process.wait(); methylation_file.close()

        ##
        ## Filter for current region's position(s)
        target_region = self.instance_params.config_dict['MethRegion'][0]
        methregion_positions = self.sequencepair_object.methylation_regions(target_region)
        np.warnings.filterwarnings('ignore') ## ignore numpy complaining about empty file if we skip header on variation_report with no results
        methylation_data = np.ndarray.tolist(np.genfromtxt(methylation_output, delimiter='\t', dtype=None, encoding='utf8', skip_header=1))
        ## HEADER:: chrom, pos, gc_content, reads, read_pp
        ## Filter all processed positions to those relevant to current methylation region
        if len(methylation_data) == 0:
            ## empty file, no variation reported for this alignment
            for position in methregion_positions:
                temp1 = ['{}{}'.format(target_region, 'NullMethylationReported'), position]; temp2 = ['0']*3
                methylation_data.append(temp1+temp2)
        else:
            ## file had variation report
            ## remove positions from source data not in the current CPG region
            ## iterate BACKWARDS to remove faster
            for i in xrange(len(methylation_data) - 1, -1, -1):
                element = methylation_data[i]
                if element[1] not in methregion_positions:
                    del methylation_data[i]

            ## check for positions missing for current CPG region data
            for position, result in itertools.izip_longest(methregion_positions, methylation_data):
                try:
                    if position != result[1]:
                        mismatch_idx = methylation_data.index(result)
                        temp1 = ['{}{}'.format(target_region, 'NullMethylationReported'), position]; temp2 = ['0']*3
                        methylation_data.insert(mismatch_idx, temp1+temp2)
                except TypeError:
                    ## variation_data had a None within iteration
                    temp1 = ['{}{}'.format(target_region, 'NullMethylationReported'), position]; temp2 = ['0']*3
                    methylation_data.append(temp1+temp2)

        ## assign to object
        if orientation == 'R1': self.sequencepair_object.set_forward_methylation(methylation_data)
        if orientation == 'R2': self.sequencepair_object.set_reverse_methylation(methylation_data)

    def determine_mapq(self, bamfi, orientation):
    	pass

    def determine_baseq(self, bamfi, orientation):
    	pass
