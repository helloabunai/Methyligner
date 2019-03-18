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
		## Indexing reference with BWA-METH
		build_subprocess = subprocess.Popen(['bwameth.py', 'index', index_copy], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		build_rawoutput = build_subprocess.communicate()
		build_stderr = build_rawoutput[1]
		build_subprocess.wait()

		##
		## Report stderr to file
		build_report = os.path.join(reference_index, 'IndexBuildReport.txt')
		report_file = open(build_report, 'w')
		report_file.write(build_stderr)
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

		forward_assembly = self.bwameth_instance(self.forward_reads, 'Aligning forward reads..','R1')
		reverse_assembly = self.bwameth_instance(self.reverse_reads, 'Aligning reverse reads..','R2')
		self.sequencepair_object.set_forward_assembly(forward_assembly)
		self.sequencepair_object.set_reverse_assembly(reverse_assembly)

	def bwameth_instance(self, inreads, feedback_string, io_index):

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
		## Output path
		log.info('{}{}{}{}'.format(clr.bold,'mth__ ',clr.end,feedback_string))
		sample_string = '{}_{}'.format(self.sample_root, io_index)
		alignment_outdir = os.path.join(self.target_output, sample_string)
		if not os.path.exists(alignment_outdir): force_mkdir(alignment_outdir)

		"""
		THREADS                     :: -t <INT>      :: CPU threads to utilise [1]
		min_seed_length             :: -k <INT>      :: minimum seed length [19]
		band_width                  :: -w <INT>      :: band width for banded alignment [100]
		seed_length_extension       :: -r <FLOAT>    :: look for internal seeds inside a seed longer than <val> [1.5]
		skip_seed_with_occurrence   :: -c <INT>      :: skip seeds with more than <val> occurrences [500]
		chain_drop                  :: -D <FLOAT>    :: drop chains shorter than <val> fraction of the overlapping chain [0.50]
		seeded_chain_drop           :: -W <INT>      :: discard chain if seeded bases shorter than <val>
		seq_match_score             :: -A <INT>      :: score for sequence match [1]
		mismatch_penalty            :: -B <INT>      :: penalty for mismatch [4]
		indel_penalty               :: -O [INT, INT] :: gap open penalites for ins/del [6,6]
		gap_extend_penalty          :: -E [INT, INT] :: penalty for extending gaps [1,1]
		prime_clipping_penalty      :: -L [INT, INT] :: 5' & 3' clipping penalty [5,5]
		unpaired_pairing_penalty    :: -U <INT>      :: penalty for unpaired read pair [17]
		"""

		## bwa-meth needs path to indvidual pre-conversion FASTA
		target = self.reference_sequence.split('/')[-1]
		vanilla_fasta = os.path.join(self.reference_sequence, '{}.fa'.format(target))
		converted_idx = os.path.join(self.reference_sequence, '{}{}'.format(target, '.fa.bwameth.c2t'))

		output_assembly = os.path.join(alignment_outdir, '{}assembly.sam'.format(io_index))
		output_file= open(output_assembly,'w')
		bwameth_process = subprocess.Popen(['bwameth.py', '--threads', str(THREADS),
			'--reference', vanilla_fasta, inreads], stdout=output_file, stderr=subprocess.PIPE)
		bwameth_stderr = bwameth_process.communicate()[1]; bwameth_process.wait()

		##
		## Generate an alignment report (i.e. console output to file)
		alignment_report = os.path.join(alignment_outdir, 'AlignmentReport.txt')
		report_file = open(alignment_report, 'w')
		report_file.write(bwameth_stderr)
		report_file.close()

		return output_assembly
