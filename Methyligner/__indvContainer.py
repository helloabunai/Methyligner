##/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

import os
import errno

class SequenceSample:
	def __init__(self):

		self.sample_label = ''
		self.instance_path = ''
		self.sample_qcpath = ''
		self.sample_alignpath = ''
		self.sample_analysispath = ''

		self.reference_idx = ''
		self.forward_fastq = ''
		self.reverse_fastq = ''

		self.initial_fastqc = None; self.initial_readcount = 0; self.initial_pbsq = None; self.initial_gcpcnt = 0
		self.postdmpx_fastqc = None; self.postdmpx_readcount = 0; self.postdmpx_pbsq = None; self.postdmpx_gcpcnt = 0
		self.posttrim_fastqc = None; self.posttrim_readcount = 0; self.posttrim_pbsq = None; self.posttrim_gcpcnt = 0

		self.forward_assembly = ''
		self.reverse_assembly = ''

		self.exception_raised = ''

	##
	## Setters
	def set_label(self, inlabel):
		self.sample_label = inlabel
	def set_instancepath(self, inpath):
		self.instance_path = inpath
	def set_qcpath(self, inpath):
		self.sample_qcpath = inpath
	def set_alignpath(self, inpath):
		self.sample_alignpath = inpath
	def set_analysispath(self, inpath):
		self.sample_analysispath = inpath

	def set_referenceidx(self, inidx):
		self.reference_idx = inidx
	def set_forwardfastq(self, inreads):
		self.forward_fastq = inreads
	def set_reversefastq(self, inreads):
		self.reverse_fastq = inreads

	def set_initial_fastqc(self, pathtofi):
		self.initial_fastqc = pathtofi
	def set_postdmpx_fastqc(self, pathtofi):
		self.postdmpx_fastqc = pathtofi
	def set_posttrim_fastqc(self, pathtofi):
		self.posttrim_fastqc = pathtofi
	def set_initial_readcount(self, count):
		self.initial_readcount = count
	def set_postdmpx_readcount(self, count):
		self.postdmpx_readcount = count
	def set_posttrim_readcount(self, count):
		self.posttrim_readcount = count
	def set_initial_pbsq(self, array):
		self.initial_pbsq = array
	def set_postdmpx_pbsq(self, array):
		self.postdmpx_pbsq = array
	def set_posttrim_pbsq(self, array):
		self.posttrim_pbsq = array
	def set_initial_gcpcnt(self, pcnt):
		self.initial_gcpcnt = pcnt
	def set_postdmpx_gcpcnt(self, pcnt):
		self.postdmpx_gcpcnt = pcnt
	def set_posttrim_gcpcnt(self, pcnt):
		self.posttrim_gcpcnt = pcnt

	def set_forward_assembly(self, pathtofi):
		self.forward_assembly = pathtofi
	def set_reverse_assembly(self, pathtofi):
		self.reverse_assembly = pathtofi

	def set_exception(self, exc):
		self.exception_raised = exc

	##
	## Getters
	def get_label(self):
		return self.sample_label
	def get_instancepath(self):
		return self.instance_path
	def get_qcpath(self):
		return self.sample_qcpath
	def get_alignpath(self):
		return self.sample_alignpath
	def get_analysispath(self):
		return self.sample_analysispath

	def get_referenceidx(self):
		return self.reference_idx
	def get_forwardfastq(self):
		return self.forward_fastq
	def get_reversefastq(self):
		return self.reverse_fastq

	def get_initial_fastqc(self):
		return self.initial_fastqc
	def get_postdmpx_fastqc(self):
		return self.postdmpx_fastqc
	def get_posttrim_fastqc(self):
		return self.posttrim_fastqc
	def get_initial_readcount(self):
		return self.initial_readcount
	def get_postdmpx_readcount(self):
		return self.postdmpx_readcount
	def get_posttrim_readcount(self):
		return self.posttrim_readcount
	def get_initial_pbsq(self):
		return self.initial_pbsq
	def get_postdmpx_pbsq(self):
		return self.postdmpx_pbsq
	def get_posttrim_pbsq(self):
		return self.posttrim_pbsq
	def get_initial_gcpcnt(self):
		return self.initial_gcpcnt
	def get_postdmpx_gcpcnt(self):
		return self.postdmpx_gcpcnt
	def get_posttrim_gcpcnt(self):
		return self.posttrim_gcpcnt

	def get_forward_assembly(self):
		return self.forward_assembly
	def get_reverse_assembly(self):
		return self.reverse_assembly

	def get_exception(self):
		return self.exception_raised

	##
	## Functions
	def generate_sampletree(self):
		for path in [self.sample_qcpath, self.sample_alignpath, self.sample_analysispath]:
			try:
				os.makedirs(path)
			except OSError as exc:
					if exc.errno == errno.EEXIST and os.path.isdir(path):pass
					else: raise
