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

class Quantification:
    def __init__(self):
    	self.sequencepair_object = sequencepair_object
    	self.instance_params = instance_params
    	self.reference_sequence = self.sequencepair_object.get_referenceidx()
    	self.forward_aln = self.sequencepair_object.get_forward_assembly()
    	self.reverse_aln = self.sequencepair_object.get_reverse_assembly()
    	self.sample_root = sequencepair_object.get_label()
    	self.target_output = sequencepair_object.get_analysispath()
