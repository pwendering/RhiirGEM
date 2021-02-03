#!/usr/bin/python3
# retrieve kinetic parameters from SABIO-RK (http://sabiork.h-its.org/)

import requests
import sys
import math
import os

def getOrgEntries(organismName):
	''' obtain all entry IDs for a given organism '''
	# URL
	ENTRYID_QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs'
	# query only with the organism name
	if organismName=="-":
		query_string="*"
	else:
		query_dict = {"Organism":'"'+organismName+'"'}
		query_string = ' AND '.join(['%s:%s' % (k,v) for k,v in query_dict.items()])
		
	query = {'format':'txt', 'q':query_string}
	# make GET request
	request = requests.get(ENTRYID_QUERY_URL, params = query)
	request.raise_for_status() # raise if 404 error
	# extract entry IDs
	entryIDs = [int(x) for x in request.text.strip().split('\n')]
	print('> %d matching entries found.' % len(entryIDs)) 
	
	return(entryIDs)
	
def getParams(entryIDs):
	''' obtain kinetic parameters for each given entry ID '''
	# URL
	PARAM_QUERY_URL = 'http://sabiork.h-its.org/entry/exportToExcelCustomizable'
	# construct request
	data_field = {'entryIDs[]': entryIDs}
	query = {'format':'tsv', 'fields[]':['ECNumber','Substrate','Organism','Parameter','Temperature','EnzymeType']}
	# make POST request
	request = requests.post(PARAM_QUERY_URL, params=query, data=data_field)
	request.raise_for_status()

	return(request.text)
	
def writeOutput(reqOutput, outFileName):
	''' write the output of the parameter request to a tab-separated file '''
	f = open(outFileName,"a+")
	f.write(reqOutput)
	f.close()

def main(organismName, outFileName):
	print("Results will be written to "+outFileName+"\n")
	if os.path.isfile(outFileName):
		res=input("Output file already exists, do you wish do overwrite it?[Y to continue]\n")
		if res is not "Y":
			exit
		
	print("Find all entry IDs for "+organismName)
	
	entryIDs = getOrgEntries(organismName)
	print("Obtain kinetic parameters\n")
	
	# separate data into bins of 1000 for faster processing
	binSize=1000
	roundVal=round(math.log(binSize,10))
	i=1
	while i*binSize <= round(len(entryIDs),-roundVal):
		binStart=(i-1)*binSize
		binEnd=i*binSize
		print("Processing IDs "+str(binStart)+"--"+str(binEnd))
		tmpEntryIDs=entryIDs[binStart:binEnd]
		responseText = getParams(tmpEntryIDs)
		writeOutput(responseText,outFileName)
		i+=1

if __name__ == "__main__":
	narg = len(sys.argv)
	
	if narg < 2:
		print("USAGE: getKineticParamsSABIORK.py organismName outFileName(opt)")
	else:
		print("\n~~~~~~~~~~~~~~~~~~~~~~~~")
		print("\tSTART")
		print("~~~~~~~~~~~~~~~~~~~~~~~~\n")
		
		organismName=sys.argv[1]
		
		if organismName == "-":
			fileSuffix = ""
		else:
			fileSuffix = "-"+"-".join(organismName.split(" "))
			
		if narg == 3:
			outFileName = sys.argv[2]
		else:
			outFileName = "kinetic-params"+fileSuffix+".tsv"
			
		main(organismName, outFileName)
