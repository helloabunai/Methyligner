##/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## Python shit
import os
import shutil
import subprocess
import logging as log

## Backend shit
from ..__backend import Colour as clr

class MethAlign:
	def __init__(self, sequencepair_object=None, instance_params=None):

		if sequencepair_object:
			print sequencepair_object

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
		## Indexing reference with bwa-build

		## TODO Methylated references require some fuckery..?

		build_subprocess = subprocess.Popen(['bwa', 'index', index_copy], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		build_rawoutput = build_subprocess.communicate()
		build_stderr = build_rawoutput[1]
		build_subprocess.wait()

		build_report = os.path.join(reference_index, 'IndexBuildReport.txt')
		report_file = open(build_report, 'w')
		report_file.write(build_stderr)
		report_file.close()

		return index_copy

	def get_index_path(self):

		return self.reference