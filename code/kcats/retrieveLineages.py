#!/usr/bin/python3
# obtain full lineages for organism names

import sys
import requests
import csv
import re
import urllib
import time

def ncbiRequest(query):
	# build the url
	query = query.split(" ")
	query = "+".join(query)
	url = "https://www.ncbi.nlm.nih.gov/taxonomy/?term="+query
	
	# try to send the request three times
	n=0
	success = False
	while not success and n<3:
		code = urllib.request.urlopen(url).getcode()
		if code == 200:
			response = urllib.request.urlopen(url).read()
			success = True
		else:
			print(str(n))
			n+=1
	pattern = re.compile("\?id=\d+")
	matches = re.search(pattern, str(response))
	# extract the taxonomy id)
	if matches:
		taxId = re.search(re.compile("\d+"), matches[0])[0]
	else:
		print("URL request not successful: "+query)
		taxId = ""
	return(taxId)

def uniprotRequest(query):
	if query == "":
		return ""
		
	baseUrl = "http://www.uniprot.org"
	dataSet = "taxonomy"
	resFormat = "tab"
	url = baseUrl+"/"+dataSet+"/?query=id:"+query+"&columns=id,lineage(all)&format="+resFormat
	attempt=0
	while attempt<10:
		attempt+=attempt
		try:
			response = requests.get(url)
			attempt=10
		except:
			print("request failed - attempts="+str(attempt))
			time.sleep(5)
	if response.ok:
		return(response.text)
	else:
		print("Something went wrong "+ str(response.status_code))
		return("")

def processUniprotOutput(result):
	try:
		return(result.split("\n")[1].split("\t")[8])
	except:
		return("")	

def main(kcatFile,outFileName):
	
	# output file with header
	fout = open(outFileName, "w")
	csvwriter = csv.writer(fout, delimiter = "\t")
	
	f = open(kcatFile, "r")
	line = f.readline()
	
	# processed names will be saved in a dictionary so they can be looked up
	# since a high level of repetition is expected
	dict = {}
	# read the file line by line
	c=0
	while line:
		line = line.strip("\n")
		line = line.split("\t")
		if c%100 == 0 and c>1:
			print("~~~ Processed "+str(c)+" lines ~~~")
			time.sleep(10) # prevents timeout from server
		# organism name
		name = line[2]
		
		if name in dict:
			res = dict[name]
		else:
			# obtain the NCBI taxonomy ID
			taxId = ncbiRequest(name)
			
			# send the request to uniprot
			res = uniprotRequest(taxId)
			
			# write the output to a new file where the third column in changed
			if res == "" or res == None:
				if name == "" or name == None:
					res = "unknown"
				else:
					res = name
			else:
				res = processUniprotOutput(res)
				
			# process the result to get the desired format
			dict[name] = res
			
		line[2] = res
		line[3] = float(line[3])
		csvwriter.writerow(line)
	
		line = f.readline()
		c+=1

	f.close()
	fout.close()
if __name__ == "__main__":
	narg = len(sys.argv)
	if  narg < 2:
		print("USAGE: retrieveLineages.py kcatFile outFileName(opt)")
	else:
		print("\n~~~~~~~~~~~~~~~~~~~~~~~~")
		print("\tSTART")
		print("~~~~~~~~~~~~~~~~~~~~~~~~\n")
		kcatFile = sys.argv[1]
	
		if narg < 3:
			outFileName = "kcat-complete-lineage.tsv"
		else:
			outFileName = sys.argv[2]
		
		main(kcatFile,outFileName)
