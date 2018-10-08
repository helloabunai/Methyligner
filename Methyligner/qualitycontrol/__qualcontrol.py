##/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## Python shit
import os
import glob
import subprocess
from fadapa import Fadapa
from multiprocessing import cpu_count
## Backend junk
from ..__backend import force_mkdir

THREADS = str(cpu_count())

class SeqQC:
	def __init__(self, sequencepair_object=None, instance_params=None, target_bin=None, stage=None):

		"""
		todo wRIte thiS haHA
		:param sequencepair_object: 
		:param instance_params: 
		:param target_bin: 
		:returns None
		"""

		self.sequencepair_object = sequencepair_object
		self.instance_params = instance_params
		self.stage = stage
		if target_bin == 'FQC': self.FastQC()
		if target_bin == 'CA': self.demultiplex()
		if target_bin == 'CATR': self.trimming()

	def FastQC(self):

		"""
		Run FastQC on target files
		Extract information from output
		Set it to object attributes as required
		:return: NoThInG
		"""

		## Target SeqQC/fastqc-stage outdir
		io_trunk = self.sequencepair_object.get_qcpath()
		target_output = os.path.join(io_trunk, self.stage)

		## Run process on specific data (init/trimmed/etc)
		fqfile = self.sequencepair_object.get_forwardfastq()
		force_mkdir(target_output)
		fastqc_process = subprocess.Popen(
			['fastqc', '--quiet', '--extract', '-t', THREADS, '-o', target_output, fqfile], stdout=subprocess.PIPE,
			stderr=subprocess.PIPE)
		fastqc_process.wait()

		## Remove ZIP of results
		for candidate in glob.glob(os.path.join(target_output,'*.zip')): os.remove(candidate)

		## Get path for fastqc_data.txt for current execution
		target_file = ''
		for root, dirs, files in os.walk(target_output):
			for name in files:
				if name.endswith('fastqc_data.txt'):
					target_file = os.path.join(root, name)

		## Number of reads present; for end-report i/o
		## Append path to FQC report so we can scrape at will
		f = Fadapa(target_file)
		stats = f.clean_data('Basic Statistics')
		pbsq = f.clean_data('Per base sequence quality')
		read_count = [x for x in stats if 'Total Sequences' in x][0][1]
		gc_pcnt = [x for x in stats if '%GC' in x][0][1]

		if self.stage == 'Initial':
			self.sequencepair_object.set_initial_readcount(read_count)
			self.sequencepair_object.set_initial_fastqc(target_file)
			self.sequencepair_object.set_initial_pbsq(pbsq)
			self.sequencepair_object.set_initial_gcpcnt(gc_pcnt)
		if self.stage == 'PostDMPX':
			self.sequencepair_object.set_postdmpx_readcount(read_count)
			self.sequencepair_object.set_postdmpx_fastqc(target_file)
			self.sequencepair_object.set_postdmpx_pbsq(pbsq)
			self.sequencepair_object.set_postdmpx_gcpcnt(gc_pcnt)
		if self.stage == 'PostTrim':
			self.sequencepair_object.set_posttrim_readcount(read_count)
			self.sequencepair_object.set_posttrim_fastqc(target_file)
			self.sequencepair_object.set_posttrim_pbsq(pbsq)
			self.sequencepair_object.set_posttrim_gcpcnt(gc_pcnt)

	def demultiplex(self):

		## Determine I/O
		dmpx_target = self.sequencepair_object.get_forwardfastq()
		file_root = dmpx_target.split('/')[-1].split('.')[0]
		target_out = os.path.join(self.sequencepair_object.get_qcpath(), 'DMPX_{}.fastq'.format(file_root))

		## Determine which argument to use with the adapter
		forward_position = self.instance_params.config_dict['demultiplex_flags']['@forward_position']
		adapter_argument = ''
		if forward_position == '5P': adapter_argument = '-g'
		if forward_position == '3P': adapter_argument = '-a'

		## Get information and execute cutadapt
		forward_adapter = self.instance_params.config_dict['demultiplex_flags']['@forward_adapter']
		error_rate = self.instance_params.config_dict['demultiplex_flags']['@error_rate']
		minimum_overlap = self.instance_params.config_dict['demultiplex_flags']['@min_overlap']
		cutadapt_subprocess = subprocess.Popen(
			['cutadapt', adapter_argument, forward_adapter, '-e', error_rate, '-O', minimum_overlap, '--discard-untrimmed',
			 '-o', target_out, dmpx_target],
			stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		cutadapt_subprocess.wait()

		## Update sequencepair's attribute for forward file to be the newly dmpx'd one
		self.sequencepair_object.set_forwardfastq(target_out)

	def trimming(self):

		## Determine I/O
		trim_target = self.sequencepair_object.get_forwardfastq()
		file_root = trim_target.split('/')[-1].split('.')[0]
		target_out = os.path.join(self.sequencepair_object.get_qcpath(), 'TRIM_{}.fastq'.format(file_root))

		## Determine arguments
		quality_threshold = self.instance_params.config_dict['trim_flags']['@quality_threshold']

		## Execute cutadapt
		cutadapt_subprocess = subprocess.Popen(
			['cutadapt', '-q', '{},{}'.format(quality_threshold,quality_threshold), '-o', target_out, trim_target],
			stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		cutadapt_subprocess.wait()

		## Update sequencepair's attribute for forward file to be the newly dmpx'd one
		self.sequencepair_object.set_forwardfastq(target_out)
