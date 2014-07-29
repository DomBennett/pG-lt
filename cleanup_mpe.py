## 29/07/2014
## D.J. Bennett
## A script for clearing up after mpe run (for development)

## Libs
import os
import shutil

def clean():
	## Remove log.txt in parent folder
	os.remove('log.txt')

	## Go through all subdirs and remove files in list
	files_to_remove = ['info.txt', 'log.txt','.paradict.p','.allrankids.p',\
	'.namesdict.p', '.genedict.p', '.terms.p']
	folders = os.listdir(os.getcwd())
	while folders:
		temp_files_to_remove = files_to_remove[:]
		while temp_files_to_remove:
			temp_file_to_remove = temp_files_to_remove.pop()
			try:
				os.remove(temp_file_to_remove)
			except OSError:
				pass

	## Go through all subdirs and remove folders in list
	folders_to_remove = ['1_names', '2_download', '3_alignment', '4_phylogeny']
	folders = os.listdir(os.getcwd())
	while folders:
		temp_folders_to_remove = folders_to_remove[:]
		while temp_files_to_remove:
			temp_folder_to_remove = temp_folders_to_remove.pop()
			try:
				shutil.rmtree(temp_folder_to_remove)
			except OSError:
				pass

if __name__ == '__main__':
	clean()