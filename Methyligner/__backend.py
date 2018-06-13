import numpy as np
import logging as log

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
			log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified data directory could not be found.'))
			trigger = True
		for fqfile in glob.glob(os.path.join(data_directory, '*')):
			if not (fqfile.endswith('.fq') or fqfile.endswith('.fastq') or fqfile.endswith('.fq.gz') or fqfile.endswith('.fastq.gz')):
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Non FastQ/GZ data detected in specified input directory.'))
				trigger = True
		reference_sequence = self.config_dict['@reference_sequence']
		if not os.path.isfile(reference_sequence):
			log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified reference sequence file could not be found.'))
			trigger = True
		if not (reference_sequence.endswith('.fa') or reference_sequence.endswith('.fas') or reference_sequence.endswith('.fasta')):
			log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified reference sequence file is not a fa/fas file.'))
			trigger = True

		##
		## Instance flag settings
		demultiplexing_flag = self.config_dict['instance_flags']['@demultiplex']
		if not (demultiplexing_flag == 'True' or demultiplexing_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Demultiplexing flag is not set to True/False.'))
			trigger = True
		sequence_qc_flag = self.config_dict['instance_flags']['@quality_control']
		if not (sequence_qc_flag == 'True' or sequence_qc_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Sequence Quality control flag is not set to True/False.'))
			trigger = True
		alignment_flag = self.config_dict['instance_flags']['@sequence_alignment']
		if not (alignment_flag == 'True' or alignment_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Sequence Alignment flag is not set to True/False.'))
			trigger = True

		##
		## Demultiplexing flag settings
		trim_adapter_base = ['A', 'G', 'C', 'T']
		if demultiplexing_flag == 'True':
			forward_adapter = self.config_dict['demultiplex_flags']['@forward_adapter']
			for charbase in forward_adapter:
				if charbase not in trim_adapter_base:
					log.error('{}{}{}{}'.format(Colour.red, 'shmthaln__d__ ', Colour.end, 'XML Config: Invalid character detected in forward_adapter demultiplexing flag.'))
					trigger = True
			forward_position = self.config_dict['demultiplex_flags']['@forward_position']
			if forward_position not in ['5P', '3P', 'AP']:
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Given demultiplexing forward adapter position invalid! [5P, 3P, AP]'))
				trigger = True

			reverse_adapter = self.config_dict['demultiplex_flags']['@reverse_adapter']
			for charbase in reverse_adapter:
				if charbase not in trim_adapter_base:
					log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Invalid character detected in reverse_adapter demultiplexing flag.'))
					trigger = True
			reverse_position = self.config_dict['demultiplex_flags']['@reverse_position']
			if reverse_position not in ['5P', '3P', 'AP']:
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Given demultiplexing reverse adapter position invalid! [5P, 3P, AP]'))
				trigger = True

			error_rate = self.config_dict['demultiplex_flags']['@error_rate']
			if not error_rate.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified error_rate is not a valid integer.'))
				trigger = True
			minimum_overlap = self.config_dict['demultiplex_flags']['@min_overlap']
			if not minimum_overlap.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified min_overlap is not a valid integer.'))
				trigger = True
			minimum_length = self.config_dict['demultiplex_flags']['@min_length']
			if not minimum_length == '':
				if not minimum_length.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified min_length is not a valid integer.'))
					trigger = True
			maximum_length = self.config_dict['demultiplex_flags']['@max_length']
			if not maximum_length == '':
				if not maximum_length.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified max_length is not a valid integer.'))
					trigger = True

		##
		## Trimming flag settings
		if sequence_qc_flag == 'True':
			trimming_type = self.config_dict['trim_flags']['@trim_type']
			if not (trimming_type == 'Quality' or trimming_type	== 'Adapter' or trimming_type == 'Both'):
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__  ', Colour.end, 'XML Config: Trimming type is not Quality/Adapter/Both.'))
				trigger = True
			quality_threshold = self.config_dict['trim_flags']['@quality_threshold']
			if not quality_threshold.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified quality threshold integer is invalid.'))
				trigger = True
			elif not int(quality_threshold) in range(0,39):
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified quality threshold integer out of range (0-38).'))
				trigger = True
			trim_adapters = ['-a','-g','-a$','-g^','-b']
			adapter_flag = self.config_dict['trim_flags']['@adapter_flag']
			if not (adapter_flag in trim_adapters):
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified trimming adapter not valid selection.'))
				trigger = True
			forward_adapter = self.config_dict['trim_flags']['@forward_adapter']
			for charbase in forward_adapter:
				if charbase not in trim_adapter_base:
					log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Invalid character detected in FW adapter sequence.'))
					trigger = True
			reverse_adapter = self.config_dict['trim_flags']['@reverse_adapter']
			for charbase in reverse_adapter:
				if charbase not in trim_adapter_base:
					log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Invalid character detected in RV adapter sequence.'))
					trigger = True
			error_tolerance = self.config_dict['trim_flags']['@error_tolerance']
			if not isinstance(float(error_tolerance), float):
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified error tolerance is not a valid float.'))
				trigger = True
			if not float(error_tolerance) in np.arange(0,1.1,0.01):
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified error tolerance is not 0.0 < x < 1.0.'))
				trigger = True

		##
		## Alignment flag settings
		if alignment_flag == 'True':
			min_seed_length = self.config_dict['alignment_flags']['@min_seed_length']
			if not min_seed_length.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified min_seed_length integer is invalid.'))
				trigger=True

			band_width = self.config_dict['alignment_flags']['@band_width']
			if not band_width.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified band_width integer is invalid.'))
				trigger=True

			seed_length_extension = self.config_dict['alignment_flags']['@seed_length_extension']
			if not isinstance(float(seed_length_extension), float):
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified seed_length_extension float is invalid.'))
				trigger=True

			skip_seed_with_occurrence = self.config_dict['alignment_flags']['@skip_seed_with_occurrence']
			if not skip_seed_with_occurrence.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified skip_seed_with_occurrence integer is invalid.'))
				trigger=True

			chain_drop = self.config_dict['alignment_flags']['@chain_drop']
			if not isinstance(float(chain_drop), float):
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified chain_drop float is invalid.'))
				trigger=True

			seeded_chain_drop = self.config_dict['alignment_flags']['@seeded_chain_drop']
			if not seeded_chain_drop.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified seeded_chain_drop integer is invalid.'))
				trigger=True

			seq_match_score = self.config_dict['alignment_flags']['@seq_match_score']
			if not seq_match_score.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified seq_match_score integer is invalid.'))
				trigger=True

			mismatch_penalty = self.config_dict['alignment_flags']['@mismatch_penalty']
			if not mismatch_penalty.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified mismatch_penalty integer is invalid.'))
				trigger=True

			indel_penalty_raw = self.config_dict['alignment_flags']['@indel_penalty']
			indel_penalty = indel_penalty_raw.split(',')
			for individual_indelpen in indel_penalty:
				if not individual_indelpen.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified indel_penalty integer(s) is(are) invalid.'))
					trigger=True

			gap_extend_penalty_raw = self.config_dict['alignment_flags']['@gap_extend_penalty']
			gap_extend_penalty = gap_extend_penalty_raw.split(',')
			for individual_gaextend in gap_extend_penalty:
				if not individual_gaextend.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified gap_extend_penalty integer(s) is(are) invalid.'))
					trigger=True

			prime_clipping_penalty_raw = self.config_dict['alignment_flags']['@prime_clipping_penalty']
			prime_clipping_penalty = prime_clipping_penalty_raw.split(',')
			for individual_prclip in prime_clipping_penalty:
				if not individual_prclip.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified prime_clipping_penalty integer(s) is(are) invalid.'))
					trigger=True

			unpaired_pairing_penalty = self.config_dict['alignment_flags']['@unpaired_pairing_penalty']
			if not unpaired_pairing_penalty.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Specified unpaired_pairing_penalty integer is invalid.'))
				trigger=True

		if trigger:
			log.error('{}{}{}{}'.format(Colour.red, 'mthaln__ ', Colour.end, 'XML Config: Failure, exiting.'))
			sys.exit(2)
		else:
			log.info('{}{}{}{}'.format(Colour.green, 'mthaln__ ', Colour.end, 'XML Config: Parsing parameters successful!'))
