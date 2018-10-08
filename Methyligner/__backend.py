##/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## Python shit
import os
import sys
import glob
import errno
import shutil
import subprocess
import numpy as np
import logging as log
from lxml import etree
from collections import defaultdict
from xml.etree import cElementTree


class Colour:

	def __init__(self):
		pass

	purple = '\033[95m'
	cyan = '\033[96m'
	darkcyan = '\033[36m'
	blue = '\033[94m'
	green = '\033[92m'
	yellow = '\033[93m'
	red = '\033[91m'
	bold = '\033[1m'
	underline = '\033[4m'
	end = '\033[0m'

class ConfigReader(object):

	"""
	The configuration file reader.
	Opens a configuration file, and if valid, converts the parameters within the file to a dictionary object,
	reader to be viewed through accessing the config_dict variable.
	"""

	def __init__(self, scriptdir, config_filename=None):

		##
		## Instance variables
		self.scriptdir = scriptdir
		self.config_filename = config_filename
		self.dtd_filename = scriptdir + "/config/config.dtd"

		##
		## Check for configuration file (just incase)
		if self.config_filename is None:
			log.error("No configuration file specified!")
		else:
			self.config_file = etree.parse(self.config_filename)

		##
		## Check config vs dtd, parse info to dictionary, validate vs ruleset
		self.validate_against_dtd()
		self.set_dictionary()
		self.validate_config()

	def validate_against_dtd(self):

		"""
		Validate input config against DTD ruleset
		i.e. confirms conformation of XML structure
		"""

		##
		## Open > etree.DTD object
		dtd_file = open(self.dtd_filename, 'r')
		dtd_object = etree.DTD(dtd_file)

		##
		## If validation fails, close the object (memory) and raise an error
		if not dtd_object.validate(self.config_file):
			dtd_file.close()
			log.error("DTD validation failure {0}: {1}".format(self.config_filename, dtd_object.error_log.filter_from_errors()[0]))
			sys.exit(2)
		dtd_file.close()

	def set_dictionary(self):

		"""
		Takes the now validated XML and extracts information from the tree into
		a python dictionary {key: value}. This dictionary will be used for variables
		within the pipeline. Recursion adapted from http://stackoverflow.com/a/9286702
		"""
		def recursive_generation(t):

			d = {t.tag: {} if t.attrib else None}
			children = list(t)

			##
			## If list was populated, create dictionary, Append keys
			if children:
				dd = defaultdict(list)

				for dc in map(recursive_generation, children):
					for k, v in dc.iteritems():
						dd[k].append(v)
				d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.iteritems()}}

			##
			## Values for key
			if t.attrib:
				d[t.tag].update(('@' + k, v) for k, v in t.attrib.iteritems())

			if t.text:
				text = t.text.strip()
				if children or t.attrib:
					if text:
						d[t.tag]['#text'] = text
				else:
					d[t.tag] = text
			return d

		##
		## Takes the formatted xml doc, puts through generator, returns dictionary
		string_repr = etree.tostring(self.config_file, pretty_print=True)
		element_tree = cElementTree.XML(string_repr)

		self.config_dict = recursive_generation(element_tree)
		self.config_dict = self.config_dict[self.config_dict.keys()[0]]

	def validate_config(self):

		"""
		Method which validates the configuration file's contents.
		If all pass, guarantees that the settings dictionary is full of valid settings!
		"""
		trigger = False

		##
		## Main configuration instance settings
		data_directory = self.config_dict['@data_dir']
		if not os.path.exists(data_directory):
			log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified data directory could not be found.'))
			trigger = True
		for fqfile in glob.glob(os.path.join(data_directory, '*')):
			if not (fqfile.endswith('.fq') or fqfile.endswith('.fastq') or fqfile.endswith('.fq.gz') or fqfile.endswith('.fastq.gz')):
				log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Non FastQ/GZ data detected in specified input directory.'))
				trigger = True
		reference_sequence = self.config_dict['@reference_sequence']
		if not os.path.isfile(reference_sequence):
			log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified reference sequence file could not be found.'))
			trigger = True
		if not (reference_sequence.endswith('.fa') or reference_sequence.endswith('.fas') or reference_sequence.endswith('.fasta')):
			log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified reference sequence file is not a fa/fasta file.'))
			trigger = True

		##
		## Instance flag settings
		demultiplexing_flag = self.config_dict['instance_flags']['@demultiplex']
		if not (demultiplexing_flag == 'True' or demultiplexing_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Demultiplexing flag is not set to True/False.'))
			trigger = True
		sequence_qc_flag = self.config_dict['instance_flags']['@quality_control']
		if not (sequence_qc_flag == 'True' or sequence_qc_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Sequence Quality control flag is not set to True/False.'))
			trigger = True
		alignment_flag = self.config_dict['instance_flags']['@sequence_alignment']
		if not (alignment_flag == 'True' or alignment_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Sequence Alignment flag is not set to True/False.'))
			trigger = True
		methylation_flag = self.config_dict['instance_flags']['@methylation_analysis']
		if not (methylation_flag == 'True' or methylation_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Methylation Analysis flag is not set to True/False'))
			trigger = True

		##
		## Demultiplexing flag settings
		trim_adapter_base = ['A', 'G', 'C', 'T']
		if demultiplexing_flag == 'True':
			forward_adapter = self.config_dict['demultiplex_flags']['@forward_adapter']
			if forward_adapter:
				for charbase in forward_adapter:
					if charbase not in trim_adapter_base:
						log.error('{}{}{}{}'.format(Colour.red, 'shmth__d__ ', Colour.end, 'XML Config: Invalid character detected in forward_adapter demultiplexing flag.'))
						trigger = True
				forward_position = self.config_dict['demultiplex_flags']['@forward_position']
				if forward_position not in ['5P', '3P', 'AP']:
					log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Given demultiplexing forward adapter position invalid! [5P, 3P, AP]'))
					trigger = True

			reverse_adapter = self.config_dict['demultiplex_flags']['@reverse_adapter']
			if reverse_adapter:
				for charbase in reverse_adapter:
					if charbase not in trim_adapter_base:
						log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Invalid character detected in reverse_adapter demultiplexing flag.'))
						trigger = True
				reverse_position = self.config_dict['demultiplex_flags']['@reverse_position']
				if reverse_position not in ['5P', '3P', 'AP']:
					log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Given demultiplexing reverse adapter position invalid! [5P, 3P, AP]'))
					trigger = True

			error_rate = self.config_dict['demultiplex_flags']['@error_rate']
			if not error_rate.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified error_rate is not a valid integer.'))
				trigger = True
			minimum_overlap = self.config_dict['demultiplex_flags']['@min_overlap']
			if not minimum_overlap.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified min_overlap is not a valid integer.'))
				trigger = True
			minimum_length = self.config_dict['demultiplex_flags']['@min_length']
			if not minimum_length == '':
				if not minimum_length.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified min_length is not a valid integer.'))
					trigger = True
			maximum_length = self.config_dict['demultiplex_flags']['@max_length']
			if not maximum_length == '':
				if not maximum_length.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified max_length is not a valid integer.'))
					trigger = True

		##
		## Trimming flag settings
		if sequence_qc_flag == 'True':
			trimming_type = self.config_dict['trim_flags']['@trim_type']
			if not (trimming_type == 'Quality' or trimming_type	== 'Adapter' or trimming_type == 'Both'):
				log.error('{}{}{}{}'.format(Colour.red, 'mth__  ', Colour.end, 'XML Config: Trimming type is not Quality/Adapter/Both.'))
				trigger = True
			## Quality
			if trimming_type == 'Quality' or trimming_type == 'Both':
				quality_threshold = self.config_dict['trim_flags']['@quality_threshold']
				if not quality_threshold.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified quality threshold integer is invalid.'))
					trigger = True
				elif not int(quality_threshold) in range(0,39):
					log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified quality threshold integer out of range (0-38).'))
					trigger = True
			
			## Adapter trimming
			trim_adapters = ['-a','-g','-a$','-g^','-b']
			adapter_flag = self.config_dict['trim_flags']['@adapter_flag']
			if trimming_type == 'Adapter' or trimming_type == 'Both':
				if not (adapter_flag in trim_adapters):
					log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified trimming adapter not valid selection.'))
					trigger = True
				forward_adapter = self.config_dict['trim_flags']['@forward_adapter']
				for charbase in forward_adapter:
					if charbase not in trim_adapter_base:
						log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Invalid character detected in FW adapter sequence.'))
						trigger = True
				reverse_adapter = self.config_dict['trim_flags']['@reverse_adapter']
				for charbase in reverse_adapter:
					if charbase not in trim_adapter_base:
						log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Invalid character detected in RV adapter sequence.'))
						trigger = True
			
			## Error tolerance
			error_tolerance = self.config_dict['trim_flags']['@error_tolerance']
			if not isinstance(float(error_tolerance), float):
				log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified error tolerance is not a valid float.'))
				trigger = True
			if not float(error_tolerance) in np.arange(0,1.1,0.01):
				log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Specified error tolerance is not 0.0 < x < 1.0.'))
				trigger = True

		##
		## Alignment flag settings
		if alignment_flag == 'True':
			placeholder = self.config_dict['alignment_flags']['@placeholder']

		##
		## Methylation Analysis flag settings
		if methylation_flag == 'True':
			placeholder = self.config_dict['analysis_flags']['@placeholder']

		if trigger:
			log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'XML Config: Failure, exiting.'))
			sys.exit(2)
		else:
			log.info('{}{}{}{}'.format(Colour.green, 'mth__ ', Colour.end, 'XML Config: Parsing parameters successful!'))

def sanitise_inputs(parsed_arguments):

	"""
	Check specified inputs exist and are correct formats
	Otherwise, trigger remains false and the inputs "pass"
	"""

	trigger = False	
	## Jobname prefix validity check
	if parsed_arguments.jobname:
		for character in parsed_arguments.jobname:
			if character is ' ' or character is '/':
				log.error('{}{}{}{}'.format(Colour.red,'mth__ ',Colour.end,'Specified Job Name has invalid characters: "', character, '"'))
				trigger = True

	## Config mode check
	if parsed_arguments.config:
		if not filesystem_exists_check(parsed_arguments.config[0]):
			log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'Specified config file could not be found.'))
			trigger = True

		for xmlfile in parsed_arguments.config:
			if not check_input_files('.xml',xmlfile):
				log.error('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'Specified config file is not an XML file.'))
				trigger = True

	return trigger

def sanitise_outputs(jobname, output_argument):

	"""
	Here we generate an output directory tree, assuming all requirements are valid
	Jobname used to be optional but i'm too lazy to give the user choice
	"""

	## Create an output directory (of specified: output + jobname) if it doesn't exist
	run_dir = ''; output_root = output_argument[0]; target_output = os.path.join(output_root, jobname)
	if not os.path.exists(target_output):
		log.info('{}{}{}{}{}'.format(Colour.bold, 'mth__ ', Colour.end, 'Creating Output with prefix: ', jobname))
		run_dir = os.path.join(output_root, jobname)
		force_mkdir(run_dir)
	## Otherwise the dir already exists..
	## make sure the user wants to overwrite before continuing
	else:
		purge_choice = ''
		while True:
			purge_choice = raw_input('{}{}{}{}'.format(Colour.bold, 'mth__ ', Colour.end,
			 'Job folder already exists. Delete existing folder? Y/N: '))
			if not (purge_choice.lower() == 'y') and not (purge_choice.lower() == 'n'):
				log.info('{}{}{}{}'.format(Colour.red, 'mth__ ', Colour.end, 'Invalid input. Please input Y or N.'))
				continue
			else:
				break

		## If the user specified to overwrite, then do so
		if purge_choice.lower() == 'y':
			log.info('{}{}{}{}{}'.format(Colour.bold, 'mth__ ', Colour.end, 'Clearing pre-existing Jobname Prefix: ', jobname))
			run_dir = os.path.join(output_root, jobname)
			if os.path.exists(run_dir):
				shutil.rmtree(run_dir, ignore_errors=True)
			force_mkdir(run_dir)
		else:
			raise Exception('User chose not to delete pre-existing Job folder. Cannot write output.')

	return run_dir

def initialise_libraries(instance_params):

	"""
	Check that all required third party binaries are on the system path!
	That way we don't need to worry about directories to stuff as a configuration option
	If they're not on the $PATH then that's the user's problem
	"""

	##
	## Subfunction for recycling code
	## Calls UNIX type for checking binaries present
	## Changed from WHICH as apparently type functions over different shells/config files
	def type_func(binary):
		binary_string = 'type {}'.format(binary)
		binary_subprocess = subprocess.Popen([binary_string], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		binary_result = binary_subprocess.communicate()
		binary_subprocess.wait()

		if 'not found' in binary_result[0] or binary_result[1]:
			log.critical('{}{}{}{}{}'.format(Colour.red,'mth__ ',Colour.end,'Missing binary: ', binary, '!'))
			raise Exception

	##
	## To determine which binaries to check for
	## AttributeError in the situation where instance_params origin differs
	## try for -c style, except AttributeError for -b style
	trigger = False
	quality_control = instance_params.config_dict['instance_flags']['@quality_control']
	alignment = instance_params.config_dict['instance_flags']['@sequence_alignment']
	analysis = instance_params.config_dict['instance_flags']['@methylation_analysis']

	if quality_control == 'True':
		try:type_func('fastqc')
		except Exception: trigger=True
		try:type_func('cutadapt')
		except Exception: trigger=True
	if alignment == 'True':
		try:type_func('bismark')
		except Exception: trigger=True
		try:type_func('bowtie2')
		except Exception: trigger=True
		try:type_func('samtools')
		except Exception: trigger=True
	if analysis == 'True':
		try:type_func('samtools')
		except Exception: trigger=True
		try:type_func('pysamstats')
		except Exception: trigger=True

	return trigger

def extract_data(input_data_directory):

	target_files = glob.glob(os.path.join(input_data_directory, '*'))
	for extract_target in target_files:
		if extract_target.lower().endswith(('.fq.gz', '.fastq.gz')):
			log.info('{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Detected compressed input data. Extracting!'))
			break

	for extract_target in target_files:
		unzipd = subprocess.Popen(['gzip', '-q', '-f', '-d', extract_target], stderr=subprocess.PIPE)
		unzipd.wait()

	return True

def sequence_pairings(data_path, instance_rundir):

	##
	## Get input files from data path
	## Sort so that ordering isn't screwy on linux
	input_files = glob.glob(os.path.join(data_path, '*'))
	sorted_input = sorted(input_files)
	sequence_pairs = []

	file_count = len(sorted_input)
	if not file_count % 2 == 0:
		log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'I/O: Non-even number of input files specified. Cannot continue without pairing!'))
		sys.exit(2)

	##
	## Optimise so code isn't recycled
	for i in range(0, len(sorted_input), 2):
		file_pair = {}
		forward_data = sorted_input[i]
		reverse_data = sorted_input[i+1]

		##
		## Check forward ends with R1
		forward_data_name = sorted_input[i].split('/')[-1].split('.')[0]
		if not forward_data_name.endswith('_R1'):
			log.error('{}{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'I/O: Forward input file does not end in _R1. ', forward_data))
			sys.exit(2)

		##
		## Check reverse ends with R2
		reverse_data_name = sorted_input[i+1].split('/')[-1].split('.')[0]
		if not reverse_data_name.endswith('_R2'):
			log.error('{}{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'I/O: Reverse input file does not end in _R2. ', reverse_data))
			sys.exit(2)

		##
		## Make Stage outputs for use in everywhere else in pipeline
		sample_root = '_'.join(forward_data_name.split('_')[:-1])
		instance_path = os.path.join(instance_rundir)
		seq_qc_path = os.path.join(instance_rundir, sample_root, 'SeqQC')
		align_path = os.path.join(instance_rundir, sample_root, 'Align')
		methyl_path = os.path.join(instance_rundir, sample_root, 'MethylationAnalysis')
		file_pair[sample_root] = [forward_data, reverse_data, instance_path, seq_qc_path, align_path, methyl_path]
		sequence_pairs.append(file_pair)

	return sequence_pairs

def filesystem_exists_check(path, raise_exception=True):

	"""
	Checks to see if the path, specified by parameter path, exists. Can be either a directory or file.
	If the path exists, True is returned. If the path does not exist, and raise_exception is set to True,
	an IOError is raised - else False is returned.
	"""

	if os.path.lexists(path):
		return True
	if raise_exception:
		log.error('{}{}{}{}'.format(Colour.red,'mth__ ',Colour.end,'Specified input path could not be found.'))
	return False

def check_input_files(input_format, input_file):

	"""
	File extension checker
	"""

	if input_file.endswith(input_format):
		return True
	return False

def force_mkdir(path):

	"""
	Some flavours of UNIX demonstrated issues with using mkdir via p2.7 os
	so fuck it, do this by force
	"""

	try: os.makedirs(path)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(path): pass
		else: raise