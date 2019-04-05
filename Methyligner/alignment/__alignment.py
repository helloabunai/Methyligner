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
		self.reference_file = reference_file
		self.target_output = target_output
		self.reference = reference_file
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
		self.reference_file = os.path.join(reference_index, self.reference.split('/')[-1])
		shutil.copy(self.reference, self.reference_file)

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

	def get_fai_faidx(self):

		return self.reference, self.reference_file

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

		forward_assembly, forward_converted = self.bwameth_instance(self.forward_reads, 'Aligning forward reads..','R1')
		reverse_assembly, reverse_converted = self.bwameth_instance(self.reverse_reads, 'Aligning reverse reads..','R2')
		self.sequencepair_object.set_forward_assembly(forward_assembly)
		self.sequencepair_object.set_reverse_assembly(reverse_assembly)
		self.sequencepair_object.set_forward_convassembly(forward_converted)
		self.sequencepair_object.set_reverse_convassembly(reverse_converted)

	def bwameth_instance(self, inreads, feedback_string, io_index):

		"""
		DoCSTRINGz
		:return:
		"""

		##
		## Tons of alignment flags
		min_seed_length = self.instance_params.config_dict['alignment_flags']['@min_seed_length']
		band_width = self.instance_params.config_dict['alignment_flags']['@band_width']
		seed_length_extension = self.instance_params.config_dict['alignment_flags']['@seed_length_extension']
		skip_seed_with_occurrence = self.instance_params.config_dict['alignment_flags']['@skip_seed_with_occurrence']
		chain_drop = self.instance_params.config_dict['alignment_flags']['@chain_drop']
		seeded_chain_drop = self.instance_params.config_dict['alignment_flags']['@seeded_chain_drop']
		seq_match_score = self.instance_params.config_dict['alignment_flags']['@seq_match_score']
		mismatch_penalty = self.instance_params.config_dict['alignment_flags']['@mismatch_penalty']
		indel_penalty = self.instance_params.config_dict['alignment_flags']['@indel_penalty']
		gap_extend_penalty = self.instance_params.config_dict['alignment_flags']['@gap_extend_penalty']
		prime_clipping_penalty = self.instance_params.config_dict['alignment_flags']['@prime_clipping_penalty']
		unpaired_pairing_penalty = self.instance_params.config_dict['alignment_flags']['@unpaired_pairing_penalty']

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

		output_assembly = os.path.join(alignment_outdir, '{}_assembly.sam'.format(sample_string))
		output_file= open(output_assembly,'w')
		bwameth_process = subprocess.Popen(['bwameth.py', '--threads', str(THREADS), '--reference', vanilla_fasta, inreads, '-k', min_seed_length,
										'-w', band_width, '-r', seed_length_extension,
										'-c', skip_seed_with_occurrence, '-D', chain_drop, '-W', seeded_chain_drop,
										'-A', seq_match_score, '-B', mismatch_penalty, '-O', indel_penalty,
										'-E', gap_extend_penalty, '-L', prime_clipping_penalty,
										'-U', unpaired_pairing_penalty], stdout=output_file, stderr=subprocess.PIPE)
		bwameth_stderr = bwameth_process.communicate()[1]; bwameth_process.wait()
		output_file.close()

		##
		## Generate an alignment report (i.e. console output to file)
		alignment_report = os.path.join(alignment_outdir, 'AlignmentReport.txt')
		report_file = open(alignment_report, 'w')
		report_file.write(bwameth_stderr)
		report_file.close()

		##
		## Convert assembly from sam to bam
		output_converted = os.path.join(alignment_outdir, '{}_converted.bam'.format(sample_string))
		converted_file = open(output_converted,'w')
		samtools_process = subprocess.Popen(['samtools', 'view' ,'-S', '-b', output_assembly],
			stdout=converted_file, stderr=subprocess.PIPE)
		samtools_stderr = samtools_process.communicate()[1]; samtools_process.wait()
		converted_file.close()

		##
		## Sort our binary assembly
		sorted_converted = os.path.join(alignment_outdir, '{}_sorted.bam'.format(sample_string))
		sorted_file = open(sorted_converted, 'w')
		sort_process = subprocess.Popen(['samtools', 'sort', output_converted],
			stdout=sorted_file, stderr=subprocess.PIPE)
		sort_stderr = sort_process.communicate()[1]; sort_process.wait()

		##
		## index the newly created sorted BAM file
		index_process = subprocess.Popen(['samtools', 'index', sorted_converted],
		 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		index_output = index_process.communicate(); index_process.wait()

		return output_assembly, sorted_converted
