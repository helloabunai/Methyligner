##/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## Python shit
import os
import shutil
import subprocess
import logging as log
import multiprocessing

## Backend shit
THREADS = multiprocessing.cpu_count()
from ..__backend import force_mkdir
from ..__backend import Colour as clr


class ReferenceIndex:
	def __init__(self, reference_file, target_output):
		self.reference = reference_file
		self.target_output = target_output
		self.reference = self.index_reference()

	def index_reference(self):

		##
		## Be paranoid, check existence/validity of reference.. again
		reference_root = self.reference.split('/')[-1].split('.')[0]
		if os.path.isfile(self.reference):
			if not (self.reference.endswith('.fa') or self.reference.endswith('.fas') or self.reference.endswith('.fasta')):
				log.critical('{}{}{}{}'.format(clr.red,'mth__ ',clr.end,'Specified reference does not exist/is not fasta.'))
		##
		## Path to store indexes for this reference
		reference_index = os.path.join(self.target_output, reference_root)
		index_copy = os.path.join(reference_index, self.reference.split('/')[-1])
		if not os.path.exists(reference_index): os.makedirs(reference_index)
		shutil.copy(self.reference, os.path.join(reference_index, self.reference.split('/')[-1]))

		##
		## Indexing reference with bowtie2-build (bismark utilises bowtie2)
		build_subprocess = subprocess.Popen(['bowtie2-build', '-f', index_copy, index_copy], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		build_rawoutput = build_subprocess.communicate()
		build_stderr = build_rawoutput[1]
		build_subprocess.wait()

		##
		## Methylation specific indexing
		meth_subprocess = subprocess.Popen(['bismark_genome_preparation', reference_index], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		meth_rawoutput = meth_subprocess.communicate()
		meth_stderr = meth_rawoutput[1]
		meth_subprocess.wait()

		##
		## Report stderr to file
		build_report = os.path.join(reference_index, 'IndexBuildReport.txt')
		report_file = open(build_report, 'w')
		report_file.write(build_stderr)
		report_file.write('\n'+meth_stderr)
		report_file.close()

		return reference_index

	def get_index_path(self):

		return self.reference

class MethAlign:
	def __init__(self, sequencepair_object=None, instance_params=None):

		self.sequencepair_object = sequencepair_object
		self.instance_params = instance_params
		self.reference_sequence = self.sequencepair_object.get_referenceidx()
		self.forward_reads = self.sequencepair_object.get_forwardfastq()
		self.reverse_reads = self.sequencepair_object.get_reversefastq()
		self.sample_root = sequencepair_object.get_label()
		self.target_output = sequencepair_object.get_alignpath()

		self.alignment_workflow()

	def alignment_workflow(self):

		forward_assembly = self.bismark_instance(self.forward_reads, 'Aligning forward reads..','R1')
		reverse_assembly = self.bismark_instance(self.reverse_reads, 'Aligning reverse reads..','R2')

	def bismark_instance(self, inreads, feedback_string, io_index):

		"""
		DoCSTRINGz
		:return: 
		"""

		##
		## Tons of alignment flags
		seed_mismatch = self.instance_params.config_dict['alignment_flags']['@seed_mismatch']
		seed_length = self.instance_params.config_dict['alignment_flags']['@seed_length']
		min_valid_insertions = self.instance_params.config_dict['alignment_flags']['@min_valid_insertions']
		max_valid_insertions = self.instance_params.config_dict['alignment_flags']['@max_valid_insertions']
		seed_extension_attempts = self.instance_params.config_dict['alignment_flags']['@seed_extension_attempts']
		reseed_attempts = self.instance_params.config_dict['alignment_flags']['@reseed_attempts']
		read_open_extend = self.instance_params.config_dict['alignment_flags']['@read_open_extend']
		rfnc_open_extend = self.instance_params.config_dict['alignment_flags']['@rfnc_open_extend']

		##
		## Bismark thread limiter
		global THREADS
		if THREADS > 4: THREADS = 4

		##
		## Output path
		log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,feedback_string))
		sample_string = '{}_{}'.format(self.sample_root, io_index)
		alignment_outdir = os.path.join(self.target_output, sample_string)
		if not os.path.exists(alignment_outdir): force_mkdir(alignment_outdir)

		"""
		THREADS                  :: -P <INT>        :: CPU threads to utilise [1]
		seed_mismatch            :: -N <INT>        :: Number of mismatches allowed in a seed during alignment [0]
		seed_length              :: -L <INT>        :: Seed substring length during alignment [20]
		min_valid_insertions     :: -I <INT>        :: Minimum insert size during alignment [0]
		max_valid_insertions     :: -X <INT>        :: Maximum insert size during alignment [500]
		seed_extension_attempts  :: -D <INT>        :: Extend seeds <int> times before failure [15]
		reseed_attempts          :: -R <INT>        :: Re-seed repetitive reads <int> times before moving on [2]
		read_open_extend         :: -rdg <INT, INT> :: Gap open/extend penalties for reads [5,3]
		rfnc_open_extend         :: -rfg <INT, INT> :: Gap open/extend penalties for reference [5,3]
		bowtie2                  :: --bowtie2       :: utilises bowtie2 aligner
		stdout                   :: n/a             :: BAM file goes to stdout
		"""

		bismark_subprocess = subprocess.Popen(['bismark', '-P', str(THREADS), '-N', seed_mismatch, '-L', seed_length,
											   '-D', seed_extension_attempts, '-R', reseed_attempts,
											   '-rdg', read_open_extend, '-rfg', rfnc_open_extend,
											   '--bowtie2', '--quiet', '-o', alignment_outdir,
											   self.reference_sequence, inreads],
											   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		bismark_rawout = bismark_subprocess.communicate(); bismark_subprocess.wait()

		##
		## Generate an alignment report (i.e. console output to file)
		alignment_report = os.path.join(alignment_outdir, 'AlignmentReport.txt')
		report_file = open(alignment_report, 'w')
		report_file.write(bismark_rawout[0])
		report_file.write(bismark_rawout[1])
		report_file.close()

		return 0




