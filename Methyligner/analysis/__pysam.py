##/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## Python shit
import os
import pysam
import shutil
import pysamstats
import subprocess
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
    	for assembly in [self.forward_aln, self.reverse_aln]:
    		self.determine_request(assembly)

    def determine_request(self, bamfi):

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

    	if quantify_variation: self.determine_variation()
    	if quantify_mapq: self.determine_mapq()
    	if quantify_baseq: self.determine_baseq()

    def determine_variation(self):
    	print 'hi'

    def determine_mapq(self):
    	print 'yo'

    def determine_baseq(self):
    	print 'howdy'
