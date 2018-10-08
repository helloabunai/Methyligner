##/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## Python shit
import os
import sys
import shutil
import argparse
import logging as log
import multiprocessing

## Backend shit
from __backend import force_mkdir
from __backend import extract_data
from __backend import ConfigReader
from __backend import Colour as clr
from __backend import sanitise_inputs
from __backend import sanitise_outputs
from __backend import sequence_pairings
from __backend import initialise_libraries
from __indvContainer import SequenceSample

## Stages
from . import qualitycontrol
from . import alignment

##GLOBALS
THREADS = multiprocessing.cpu_count()

class Methyligner:
	def __init__(self):
		
		"""
		TODO DOCSTRINGS
		"""

		##
		## Argument parser from CLI
		self.parser = argparse.ArgumentParser(prog='methyligner', description='Methyligner: Alignment/Analysis pipeline for Methylated DNA sequences.')
		self.parser.add_argument('-v', '--verbose', help='Verbose output mode. Setting this flag enables verbose output. Default: off.', action='store_true')
		self.parser.add_argument('-c', '--config', help='Pipeline config. Specify a directory to your ArgumentConfig.xml file.', nargs=1, required=True)
		self.parser.add_argument('-t', '--threads', help='Thread utilisation. Typically only alters third party alignment performance. Default: system max.', type=int, choices=xrange(1, THREADS+1), default=THREADS)
		self.parser.add_argument('-j', '--jobname', help='Customised folder output name. If not specified, defaults to normal output naming schema.', type=str, required=True)
		self.parser.add_argument('-o', '--output', help='Output path. Specify a directory you wish output to be directed towards.', metavar='output', nargs=1, required=True)
		self.args = self.parser.parse_args()
		self.header = ''

		##
		## Set verbosity for CLI output
		if self.args.verbose:
			log.basicConfig(format='%(message)s', level=log.DEBUG)
			log.info('{}{}{}{}'.format(clr.bold, 'mth__ ', clr.end, 'Methyligner: Alignment/Analysis pipeline for Methylated DNA sequences.'))
			log.info('{}{}{}{}'.format(clr.bold, 'mth__ ', clr.end, 'alastair.maxwell@glasgow.ac.uk\n'))
		else:
			log.basicConfig(format='%(message)s')

		##
		## Check we're on python 2.7.13 or greater (dependencies require)
		if not (sys.version_info[0] == 2 and sys.version_info[1] == 7):
			if sys.version_info[2] < 13:
				current_user_version = '{}.{}.{}'.format(sys.version_info[0], sys.version_info[1], sys.version_info[2])
				log.error('{}{}{}{}{}.'.format(clr.red, 'mth__ ', clr.end, 'Methyligner requires python2 2.7.13 or later!'
				' You are using: ', current_user_version))
				sys.exit(2)

		##
		## Check inputs are valid, func returns True if errors detected
		if sanitise_inputs(self.args):
			log.error('{}{}{}{}'.format(clr.red, 'mth__ ', clr.end, 'Error with specified input(s). Exiting!'))
			sys.exit(2)
		##
		## Attempt generation of output dirtree
		try: self.instance_rundir = sanitise_outputs(self.args.jobname, self.args.output)
		except Exception, e: log.error('{}{}{}{}'.format(clr.red, 'mth__ ', clr.end, e)); sys.exit(2)

		##
		## Set up config dictionary of input parameters
		## Copy configuration file to instance output folder (for reproducability of runs)
		self.configfile = self.args.config[0]; script_path = os.path.dirname(__file__)
		instance_configuration = os.path.join(self.instance_rundir, 'UtilisedConfiguration.xml')
		shutil.copyfile(self.configfile, instance_configuration)
		self.instance_params = ConfigReader(script_path, self.configfile)
		self.instance_params.config_dict['JobName'] = self.args.jobname

		##
		## Check system $PATH for required third party binaries
		if initialise_libraries(self.instance_params):
			log.error('{}{}{}{}'.format(clr.red, 'mth__ ', clr.end, 'Detected missing library from system/$PATH. Exiting!'))
			sys.exit(2)
		else:
			log.info('{}{}{}{}'.format(clr.green,'mth__ ', clr.end, 'Required libraries present, assuming OK!\n'))

		##
		## Set up instance-wide applicable files
		## TODO reset when applicable to certain stages
		self.index_path = ''; self.reference_indexes = [];
		self.instance_results = ''; self.instance_graphs = ''

		##
		## Workflow time!
		self.sequence_workflow()

		##
		## Finished!
		log.info('{}{}{}{}'.format(clr.green, 'mth__ ', clr.end, 'Methyligner pipeline completed; exiting!'))

	def sequence_workflow(self):
		"""
		Workflow for when Methyligner is being ran in config mode..
		Behaviours are tailored based on information extracted from the specified config XML file
		:return: None
		"""		

		##
		## Data location
		instance_inputdata = self.instance_params.config_dict['@data_dir']

		##
		## Pre-processing; check for compressed data
		if not extract_data(instance_inputdata):
			log.error('{}{}{}{}'.format(clr.red, 'mth__ ', clr.end, 'Error during file extraction. Please check data!'))
			sys.exit(2)

		##
		## Pre-processing; index references
		if self.instance_params.config_dict['instance_flags']['@sequence_alignment']:
			log.info('{}{}{}{}'.format(clr.bold,'mth__ ',clr.end,'Indexing reference(s) before initialising sample pair cycle..'))
			self.index_path = os.path.join(self.instance_rundir,'Indexes'); force_mkdir(self.index_path)
			reference_sequence = self.instance_params.config_dict['@reference_sequence']
			reference_indexes = alignment.ReferenceIndex(reference_sequence, self.index_path).get_index_path()
			self.reference_indexes = reference_indexes

		data_pairs = sequence_pairings(instance_inputdata, self.instance_rundir)
		for i in range(len(data_pairs)):
			for seqpair_lbl, seqpair_dat in data_pairs[i].iteritems():
				################################################
				## Pre stage! Sample object/Tree generation.. ##
				################################################
				log.info('\n{}{}{}{}{}/{} ({})'.format(clr.bold, 'mth__ ', clr.end, 'Processing sequence pair: ',
													 str(i + 1), str(len(data_pairs)), seqpair_lbl))
				current_seqpair = SequenceSample()
				current_seqpair.set_label(seqpair_lbl)
				current_seqpair.set_instancepath(seqpair_dat[2])
				current_seqpair.set_qcpath(seqpair_dat[3])
				current_seqpair.set_alignpath(seqpair_dat[4])
				current_seqpair.set_analysispath(seqpair_dat[5])
				current_seqpair.set_referenceidx(self.reference_indexes)
				current_seqpair.set_forwardfastq(seqpair_dat[0])
				current_seqpair.set_reversefastq(seqpair_dat[1])
				current_seqpair.generate_sampletree()

				#########################################
				## Stage 1 (A) FastQC on initial data. ##
				#########################################
				log.info('{}{}{}{}'.format(clr.bold, 'mth__ ', clr.end, 'Determining input read quality..'))
				try:
					qualitycontrol.SeqQC(current_seqpair, self.instance_params, target_bin='FQC', stage='Initial')
					log.info('{}{}{}{}'.format(clr.green, 'mth__ ', clr.end, 'Complete!'))
				except Exception, e:
					current_seqpair.set_exception('SeqQC-FQCInit')
					log.info('{}{}{}{}{}: {}\n'.format(clr.red,'mth__ ',clr.end,'SeqQC failure on ',seqpair_lbl,str(e)))
					continue

				##################################################
				## Stage 1 (B) Demultiplex via forward primer.. ##
				##################################################
				log.info('{}{}{}{}'.format(clr.bold, 'mth__ ', clr.end, 'Demultiplexing sequence reads..'))
				try:
					qualitycontrol.SeqQC(current_seqpair, self.instance_params, target_bin='CA')
					log.info('{}{}{}{}'.format(clr.green, 'mth__ ', clr.end, 'Complete!'))
				except Exception, e:
					current_seqpair.set_exception('SeqQC-DMPX')
					log.info('{}{}{}{}{}: {}\n'.format(clr.red,'mth__ ',clr.end,'SeqQC failure on ',seqpair_lbl,str(e)))
					continue

				####################################################
				## Stage 1 (C) FastQC on DMPX to determine CpG5.. ##
				####################################################
				log.info('{}{}{}{}'.format(clr.bold, 'mth__ ', clr.end, 'Processing demultiplexed reads for quality..'))
				try:
					qualitycontrol.SeqQC(current_seqpair, self.instance_params, target_bin='FQC', stage='PostDMPX')
					log.info('{}{}{}{}'.format(clr.green, 'mth__ ', clr.end, 'Complete!'))
				except Exception, e:
					current_seqpair.set_exception('SeqQC-DMPX')
					log.info('{}{}{}{}{}: {}\n'.format(clr.red,'mth__ ',clr.end,'SeqQC failure on ',seqpair_lbl,str(e)))
					continue

				####################################
				## Stage 1 (D) Quality trimming.. ##
				####################################
				log.info('{}{}{}{}'.format(clr.bold, 'mth__ ', clr.end, 'Trimming by sequence PHRED score..'))
				try:
					qualitycontrol.SeqQC(current_seqpair, self.instance_params, target_bin='CATR')
					log.info('{}{}{}{}'.format(clr.green, 'mth__ ', clr.end, 'Complete!'))
				except Exception, e:
					current_seqpair.set_exception('SeqQC-DMPX')
					log.info('{}{}{}{}{}: {}\n'.format(clr.red,'mth__ ',clr.end,'SeqQC failure on ',seqpair_lbl,str(e)))
					continue

				###########################################
				## Stage 1 (E) FastQC on post-PHRED trim ##
				###########################################
				log.info('{}{}{}{}'.format(clr.bold, 'mth__ ', clr.end, 'Processing PHRED filtered reads for quality..'))
				try:
					qualitycontrol.SeqQC(current_seqpair, self.instance_params, target_bin='FQC', stage='PostTrim')
					log.info('{}{}{}{}'.format(clr.green, 'mth__ ', clr.end, 'Complete!'))
				except Exception, e:
					current_seqpair.set_exception('SeqQC-DMPX')
					log.info('{}{}{}{}{}: {}\n'.format(clr.red,'mth__ ',clr.end,'SeqQC failure on ',seqpair_lbl,str(e)))
					continue

				print '\n>>HELLO'
				print 'Init: ', current_seqpair.get_initial_readcount(), current_seqpair.get_initial_gcpcnt()
				print 'PostDMPX: ', current_seqpair.get_postdmpx_readcount(), current_seqpair.get_postdmpx_gcpcnt()
				print 'PostTrim: ', current_seqpair.get_posttrim_readcount(), current_seqpair.get_posttrim_gcpcnt()

				###############################################
				## Stage 2 (lol finally) Bismark alignment!! ##
				###############################################

def main():
	try:
		Methyligner()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))
		sys.exit(2)