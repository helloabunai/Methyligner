##/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

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

	##
	## Functions
	def generate_sampletree(self):
		for path in [self.sample_qcpath, self.sample_alignpath, self.sample_analysispath]:
			try:
				os.makedirs(path)
			except OSError as exc:
					if exc.errno == errno.EEXIST and os.path.isdir(path):pass
					else: raise