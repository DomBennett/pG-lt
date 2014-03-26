#!/usr/bin/python
## MPE Entrez tools
## D.J. Bennett
## 24/03/2014

## Packages
import re, time, sys, os, random
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

## Globals
max_check = 4
download_counter = 0

## Functions
def eSearch(term, retStart=0, retMax=1, usehistory="n", db = "nucleotide"):
	"""Use Entrez.esearch to search a term in an NCBI database.

	Arguments:
	 term = string of term used in search
	 retStart = minimum returned ID of matching sequences IDs
	 retMax = maximum returned ID of matching sequences IDs
	 usehistory = record search in NCBI database ("y" or "n")
	 db = NCBI database
	 
	Return:
	 dictionary"""
	finished = 0
	global download_counter
	while finished <= max_check:
		if download_counter > 1000:
			print "Download counter hit. Waiting for 60 seconds ..."
			download_counter = 0
			time.sleep(60)
		try:
			if db is "nucleotide":
				handle = Entrez.esearch(db="nucleotide",term=term, usehistory=usehistory,\
								retStart=retStart, retMax=retMax, retmode="text")
				results = Entrez.read(handle)
				handle.close()
			elif db is "taxonomy":
				handle = Entrez.esearch(db = "taxonomy", retStart = retStart, retmax = retMax,\
					term = term)
				results = Entrez.read(handle)
				handle.close()
			else:
				print "Invalid db argument!"
				break
			download_counter += 1
			finished = max_check + 1
		except:
			if finished == 0:
				print "!!!Server error checking", term, " - retrying..."
				time.sleep(10)
			elif finished == max_check:
				print "!!!!!!Unreachable. Returning nothing."
				return()
			else:
				finished += 1
				time.sleep(10)
	return results

def eFetch(ncbi_id, db = "nucleotide"):
	"""Download NCBI record(s) using ID number(s).

	Arguments:
	 ncbi_id = sequence identifier (list or string)

	Return:
	 SeqRecord object"""
	finished = 0
	global download_counter
	while finished <= max_check:
		if download_counter > 1000:
			print "Download counter hit. Waiting for 60 seconds ..."
			download_counter = 0
			time.sleep(60)
		try:
			if db is "nucleotide":
				handle = Entrez.efetch(db = "nucleotide", rettype = 'gb',\
						       retmode="text", id = ncbi_id)
				results_iter = SeqIO.parse(handle, 'gb')
				results = [x for x in results_iter]
				handle.close()
			elif db is "taxonomy":
				handle = Entrez.efetch(db = "taxonomy", id = ncbi_id, retmode = "xml")
				results = Entrez.read(handle)
				handle.close()
			else:
				print "Invalid db argument!"
				break
			download_counter += len(ncbi_id)
			finished = max_check + 1
		except ValueError: # if parsing fails, value error raised
			handle.close()
			results = ()
			finished = max_check + 1
		except:
			if finished == 0:
				print "!!!Server error - retrying..."
				finished += 1
				time.sleep(10)
			if finished == max_check:
				print "!!!!!!Unreachable. Returning nothing."
				return ()
	return results

def findChildren(taxid, target = 100, next = False):
	"""
	Return all decendant genera (or below) of a taxonmic ID.

	Args:
	 taxid = taxonomic ID
	 target = the target number of children returned (default 100)
	 next = stop at all children in the rank below given id's rank

	Returns:
	 taxid
	"""
	def findNext(frecord):
		term = "{0}[Next Level]".format(frecord[0]['ScientificName'])
		count = eSearch(term, db = "taxonomy")["Count"]
		srecord = eSearch(term, db = "taxonomy", retMax = count)
		return srecord['IdList']
	def findTillTarget(taxids):
		res = []
		taxids = random.sample(taxids, len(taxids))
		while len(taxids) > 0:
			if len(res) > target:
				break
			taxid = taxids.pop()
			frecord = eFetch(taxid, db = "taxonomy")
			print frecord[0]['Rank']
			if frecord[0]['Rank'] in target_ranks:
				res.append(taxid)
			else:
				res.extend(findTillTarget(findNext(frecord)))
		return res
	taxid = str(taxid)
	target_ranks = ['genus', 'subgenus', 'species', 'subspecies']
	if next:
		frecord = eFetch(taxid, db = "taxonomy")
		return findNext(frecord)
	else:
		return findTillTarget([taxid])

if __name__ == '__main__':
	pass