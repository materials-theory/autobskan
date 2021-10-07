# coding: utf-8

__author__ = "Giyeok Lee"
__email__ = "lgy4230@yonsei.ac.kr"
__date__ = "Oct 08, 2021"
__maintainer__ = "Giyeok Lee"
__version__ = "1.0.0"
__copyright__ = "Copyright (c) Materials Theory Group @ Yonsei University (2021)"


import numpy as np
import os
import glob
from ase.io.vasp import read_vasp, write_vasp
from tqdm import tqdm

import autobskan.input.input
from autobskan.image import stmplot, post_processing, AR


def main():
	if os.path.exists("bskan.in"):
		bskan_input = autobskan.input.input.Bskan_input("bskan.in")
	else:
		# Using argparse can be another options.. But it could be somewhat messy in that case
		raise IOError("plz make bskan.in file for CLI, or autobskan-gui command for GUI")

	try:
		str_vasp = read_vasp(bskan_input.poscar)
		bskan_input.poscar = str_vasp
		bskan_input.gamma = str_vasp.cell.cellpar()[-1]
	except:
		# raise IOError("wrong input structure")
		print("No input structures. Gamma will be selected as you set. Default=90")

	if bskan_input.atom_addinfo != None:
		with open(bskan_input.atom_addinfo, 'r') as addinfo:
			bskan_input.atom_addinfo = addinfo.readlines()

	if bskan_input.mode == "IMAGE":
		if isinstance(bskan_input.current, str):
			if bskan_input.current == "ALL" :
				bskan_input.current = "*"
			current_files = glob.glob(bskan_input.current)
		else:
			current_files = bskan_input.current

		for current_file in tqdm(current_files, desc=f"Total {len(current_files)} files"):
			# try:
			current = stmplot.Current(current_file)
			if not os.path.exists(current_file+"_images"):
				os.makedirs(current_file+"_images")
			os.chdir(current_file+"_images")

			# Automatic selection of ISOSURFACE values
			if bskan_input.iso_auto == "LOGSCALE":
				iso_max, iso_min = np.floor(np.log10(current.iso_max)), np.ceil(np.log10(current.iso_min))
				bskan_input.iso = 10 ** np.arange(iso_min, iso_max+1)
			elif bskan_input.iso_auto:
				try:
					if not isinstance(bskan_input.iso, list):
						bskan_input.iso = list(np.linspace(current.iso_max, current.iso_min, int(bskan_input.iso)))
				except:
					raise IOError("wrong input of ISO in bskan.in")
			else: # ISO_AUTO = FALSE
				pass

			if not os.path.exists("raw_images"):
				os.makedirs("raw_images")
			stmplot.main(current, bskan_input, "raw_images", plot_atoms=True)
			if bskan_input.post_processing:
				post_processing.main("raw_images", bskan_input)
			os.chdir("../")
			# except:
			# 	print(f'{current_file} is not a regular bSKAN output file.')
			# 	continue

	elif bskan_input.mode=="POST_PROCESSING":
		post_processing.main() #--> gyb랑 동일하게 사용할 수 있도록!

# elif bskan_input.mode == "CALCULATION":
# else:
## TEST mode

if __name__ == "__main__":
	main()
