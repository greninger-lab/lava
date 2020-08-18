# install script for all dependencies except for annovar
import platform
import sys
import subprocess
import os
import argparse
import distutils.spawn

if __name__ == '__main__':
	# TODO: have a flag for install and a flag for checking dependencies in the same way that lava does 
	parser = argparse.ArgumentParser(description='LAVA install script. Once you have Python and an internet connection '
		'and Brew if running on a Mac, use this script to install almost all dependencies')
	parser.add_argument('-i',action='store_true', help='Use -i to download and install all dependencies')
	parser.add_argument('-c',action='store_true', help ='Use -c to check that all dependencies are available to LAVA')
	gatk_version = 'gatk-4.0.11.0/'

	args = parser.parse_args()

	if args.i == False and args.c == False:
		print('Try running python install.py -h for some help')
		sys.exit(1)

	if args.i:
		print('Please note that this script assumes you have python installed on your system and that you are executing this '
			'file from the main lava folder.')
		print('If running on Mac OS you also need to have brew properly installed.')
		print
		print('Currently you are running this install script from:')
		subprocess.call('pwd', shell=True)
		print

		# For super fresh OS systems
		if platform.system() == 'Darwin':
			print('Installing wget...')
			subprocess.call('brew list wget || brew install wget', shell = True)

		# Either install pip or upgrade it to the latest version 
		print('Installing pip...')
		if not os.path.isfile('get-pip.py'):
			subprocess.call('wget https://bootstrap.pypa.io/get-pip.py', shell=True)
			subprocess.call('python get-pip.py --user', shell=True)

		# Installing python modules
		print('Installing python modules...')
		print('First we\'re going to update pip and some other setuptools.')

		print('Installing biopython...')
		subprocess.call('python3 -m pip install biopython --user', shell=True)
		print('Installing numpy...')
		subprocess.call('python3 -m pip install numpy --user ', shell=True)
		print('Installing pandas...')
		subprocess.call('python3 -m pip install pandas --user', shell=True)
		print('Installing bokeh...')
		subprocess.call('python3 -m pip install bokeh --user', shell=True)
		print
		print('Python modules installed!')
		print

		
		# TODO: add detection to only do this if picard.jar isn't found in the main lava folder 
		if not os.path.isfile('picard.jar'):
			print('Downloading Picard...')
			subprocess.call('wget https://github.com/broadinstitute/picard/releases/download/2.18.15/picard.jar', shell=True)
		if not os.path.isdir('gatk-4.0.11.0'):
			print('Downloading GATK and installing GATK...')
			subprocess.call('wget https://github.com/broadinstitute/gatk/releases/download/4.0.11.0/gatk-4.0.11.0.zip', shell=True)
			subprocess.call('unzip gatk-4.0.11.0.zip', shell=True)
			subprocess.call('rm gatk-4.0.11.0.zip', shell=True)
		if not os.path.isfile('VarScan'):
			print('Downloading VarScan...')
			subprocess.call('wget --no-check-certificate https://sourceforge.net/projects/varscan/files/latest/download', shell=True)
			subprocess.call('mv download VarScan', shell=True)

		if platform.system() == 'Linux':
			if not os.path.isfile('gff3ToGenePred'):
				subprocess.call('wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred', shell=True)
			print
			print("Some tools require root access.")
			subprocess.call('sudo apt-get install bedtools samtools bwa mafft bcftools', shell=True)
		elif platform.system() == 'Darwin':
			if not os.path.isfile('gff3ToGenePred'):
				subprocess.call('wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/gff3ToGenePred', shell = True)
			subprocess.call('brew list bedtools || brew install bedtools', shell=True)
			subprocess.call('brew list samtools || brew install samtools', shell=True)
			subprocess.call('brew list bwa || brew install bwa', shell=True)
			subprocess.call('brew list mafft || brew install mafft', shell= True)
			subprocess.call('brew list bcftools || brew install bcftools', shell=True)
		subprocess.call('chmod +x gff3ToGenePred', shell =True)

	# This way we support doing both -i and -c one after the other to install and check in one execution 
	if args.c:
		# if this is still 0 at the end of this block then we know there's nothing fishy about the environment 
		error_code = 0 
		annovar_error = 0
		try:
			from Bio.Seq import Seq
		except ImportError:
			print('Biopython not correctly installed! - Specifically Bio.Seq is not available to Python')
			error_code += 1

		try: 
			from Bio.Blast import NCBIWWW 
		except ImportError:
			print('Biopython not correctly installed! - Specifically Bio.Blast is not available')
			error_code += 1

		try:
			from Bio import Entrez
		except ImportError:
			print('Biopython not correctly installed! - Specifically Entrez module')
			error_code += 1

		try: 
			import numpy
		except ImportError:
			print('Numpy python module not correctly installed!')
			error_code += 1

		try:
			import pandas
		except ImportError:
			print('Pandas python module not correctly installed!')
			error_code += 1
			
		try:
			import bokeh
		except ImportError:
			print('bokeh python module not correctly installed!')
			error_code += 1


		if not os.path.isfile('./picard.jar'):
			print('LAVA cannot find picard.jar - check that you are executing this script from the main LAVA directory')
			error_code += 1

		if not os.path.isdir('./' + gatk_version):
			print('LAVA cannot find ' + gatk_version + ' - this error can be caused by not executing this script from the main LAVA '
				'directory. This error can also be caused by downloading a different version of GATK. Check your main LAVA folder and '
				'if you see a folder that says something that starts with \'gatk-\' but is NOT \'' + gatk_version + '\' then you need to'
				' either go and download ' + gatk_version + ' or rename whatever gatk you have in your folder to ' + gatk_version)
			error_code += 1

		if not os.path.isfile('./VarScan'):
			print('LAVA cannot find VarScan in the main LAVA folder.')
			error_code += 1

		if not os.path.isfile('./gff3ToGenePred'):
			print('gff3ToGenePred not in the main LAVA folder. NOTE: if you are receiving errors related to this file but you have '
				'verified that it is in the main LAVA folder you might try running \'chmod +x gff3ToGenePred\' from inside the main LAVA '
				'directory')
			error_code += 1

		if not os.path.isfile('./retrieve_seq_from_fasta.pl'):
			print('retrieve_seq_from_fasta.pl not available in the main LAVA directory - this means that you have yet to download ANNOVAR '
				'and place it into the main LAVA directory. To fix this error check out the readme or simply go to '
				'http://www.openbioinformatics.org/annovar/annovar_download_form.php and request acsess, download, unzip and place all '
				'files into the main LAVA directory')
			error_code += 1
			annovar_error =1
		if not os.path.isfile('./convert2annovar.pl'):
			print('Another ANNOVAR scirpt (convert2annovar.pl)is missing, to fix download and unzip all ANNOVAR files to main LAVA directory')
			error_code += 1
			annovar_error = 1
		if not os.path.isfile('annotate_variation.pl'):
			print('Another ANNOVAR script (annotate_variation.pl) is missing to fix download and unzip all ANNOVAR files to main LAVA directory')
			error_code += 1
			annovar_error = 1

		bwa_run = distutils.spawn.find_executable('bwa')
		if bwa_run == None:
			print('bwa not properly installed, to fix this either run brew install bwa (Mac) or apt-get install bwa (Linux)')
			error_code += 1

		sam_run = distutils.spawn.find_executable('samtools')
		if sam_run == None:
			print('samtools not properly installed! To fix either run brew install samtools (Mac) or apt-get install samtools (Linux)')
			error_code += 1

		bcf_run = distutils.spawn.find_executable('bcftools')
		if bcf_run == None:
			print('bcftools not properly installed! To fix either run brew install bcftools (Mac) or apt-get install bcftools (Linux)')
			error_code += 1

		bed_run = distutils.spawn.find_executable('bedtools')
		if bed_run == None:
			print('bedtools not properly installed! To fix either run brew install bedtools (Mac) or apt-get install bedtools (Linux)')
			error_code += 1

		maft_run = distutils.spawn.find_executable('mafft')
		if maft_run == None:
			print('mafft not properly installed! to fix either run brew install mafft (mac) or apt-get install mafft (linux)')
			error_code += 1

		java_run = distutils.spawn.find_executable('java')
		if java_run == None:
			print('Something is weird with your java installation - make sure you have a JRE that can run -jar files properly '
				'installed check out https://docs.oracle.com/javase/8/docs/technotes/guides/install/install_overview.html')
			error_code += 1
		if error_code == 3 and annovar_error ==1:
			print('Everything working except ANNOVAR which is normal for first time installation - check out the error messeges or the readme for how to install ANNOVAR')
		if error_code == 0:
			print
			print('All dependencies working properly! Time to do some longitudinal analysis of viral alleles! :DDDDDD')
			print
