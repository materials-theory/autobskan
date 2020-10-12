import argparse
import shutil
import glob
import os
import matplotlib.pyplot as plt
import numpy as np

def dirlist():
    alls = glob.glob("*")
    direc = [x for x in alls if os.path.isdir(x)]
    return direc

############# CONDUCT (engine) #############
# 1. just a simple copy
def conduct(wanted, savedir, savename):
	shutil.copy(wanted, savedir + savename)

#############################################

def algo(wanted, savedir, start_dir, separator="_"):
	for i in dirlist():
		os.chdir(i)
		if wanted in glob.glob("*"):
			path = os.getcwd()
			cor_path = path[len(start_dir)+1:]
			cor_path = separator.join(cor_path.split("/"))
			conduct(wanted, savedir, cor_path)
		if len(dirlist()) != 0:
			algo(wanted, savedir, start_dir, separator)
		os.chdir("../")

def main():
	pars = argparse.ArgumentParser()
	pars.add_argument('-t', type = str, default='CURRENT', help='target file. conduct function will be performed only when target file exists')
	pars.add_argument('-o', type = str, default='collected_current', help='file will be collected in this directory (location : current location)')
	pars.add_argument('-sep', type = str, default='_', help='the paths of each target files.')
	args = pars.parse_args()
	wanted, savedir_name, sep = args.t, args.o, args.sep
	savedir = os.getcwd()+"/"+savedir_name+"/"
	os.makedirs(savedir)
	start_dir = os.getcwd()

	# wanted files which is in subdirectories
	algo(wanted, savedir, start_dir, sep)

	# wanted files which is in current directories
	if wanted in glob.glob("*"):
		path = os.getcwd()
		cor_path = path[len(start_dir)+1:]
		cor_path = sep.join(cor_path.split("/"))
		conduct(wanted, savedir, cor_path)


if __name__=="__main__":
	main()
