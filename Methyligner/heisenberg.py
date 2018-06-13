import os
import multiprocessing

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
		self.parser.add_argument('-j', '--jobname', help='Customised folder output name. If not specified, defaults to normal output naming schema.', type=str)
		self.parser.add_argument('-o', '--output', help='Output path. Specify a directory you wish output to be directed towards.', metavar='output', nargs=1, required=True)
		self.args = self.parser.parse_args()
		self.header = ''

		##
		## Set verbosity for CLI output
		if self.args.verbose:
			log.basicConfig(format='%(message)s', level=log.DEBUG)
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Methyligner: Alignment/Analysis pipeline for Methylated DNA sequences.'))
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'alastair.maxwell@glasgow.ac.uk\n'))
		else:
			log.basicConfig(format='%(message)s')

def main():
	try:
		Methyligner()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))

		sys.exit(2)