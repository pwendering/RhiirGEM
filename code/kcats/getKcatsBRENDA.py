#!/usr/bin/python3
# retrieve kcat values from BRENDA
import sys
from zeep import Client
import hashlib
import csv
import time

def processBRENDAoutput(entry):
	line = [];
	line.append(entry["ecNumber"])
	line.append(entry["substrate"])
	line.append(entry["organism"])
	line.append(entry["turnoverNumber"])
	line.append(entry["commentary"])
	return(line)
		
	
def append2KcatFile(csvwriter, result):
	for entry in result:
		line = processBRENDAoutput(entry)
		csvwriter.writerow(line)	

def getECWithTurnoverNumber(client,username,password):
	parameters = (username,password)
	ecNumbers=client.service.getEcNumbersFromTurnoverNumber(*parameters)
	return(ecNumbers)


def main(organismName, outFileName, username, password):
	
	if organismName=="-":
		organismName=""

	# try three times to get the client
	success = False
	n=0
	while n<3 and not success:
		try:
			wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
			client = Client(wsdl)
			success = True
			print("Client set up ("+str(n+1)+" attempt(s)).")
		except:
			n+=1
			
	ecNumbers = getECWithTurnoverNumber(client,username,password)

	# open output file
	print("write output to "+outFileName+"\n")
	fout = open(outFileName, "w", newline = '')
	writer = csv.writer(fout, delimiter = "\t")
	
	# write the head line
	writer.writerow(["ec_number","substrate","organism","kcat","commentary"])

	# loop over obtained E.C. numbers
	print("sending requests for "+str(len(ecNumbers))+" E.C. numbers")
	print("first entry: "+ecNumbers[0])
	print("last entry: "+ecNumbers[len(ecNumbers)-1])
	
	for i in range(0,len(ecNumbers)):
		
		success = False
		n=0
		while n<10 and not success:
			# send the request
			try:
				parameters = (username,password,"ecNumber*"+ecNumbers[i],
					"organism*"+organismName,"turnoverNumberMaximum*", "substrate*", "commentary*",
					"turnoverNumber*", "ligandStructureId*", "literature*")
				
				resultString = client.service.getTurnoverNumber(*parameters)
				# write the output to file
				append2KcatFile(writer, resultString)
				if resultString is not "":
					success = True
			except:
				time.sleep(5)
				n+=1
				print("Request attempt "+str(n)+" failed")
		if i % 10 == 0 and i > 0:
			time.sleep(0.5)
		if i % 100 == 0 and i > 0:
			print("done with "+str(i))
			time.sleep(5)

	fout.close()

if __name__ == "__main__":
	narg = len(sys.argv)
	if narg < 2:
		print("USAGE: getKcatsBRENDA.py organismName outFileName(opt) username(opt) password(opt)")
	else:
		print("\n~~~~~~~~~~~~~~~~~~~~~~~~")
		print("\tSTART")
		print("~~~~~~~~~~~~~~~~~~~~~~~~\n")
		
		organismName = sys.argv[1]
		
		if organismName=="-":
			fileSuffix=""
		else:
			fileSuffix="-"+"-".join(organismName.split(" "))
		
		if narg < 3:
			outFileName = "kcats-BRENDA"+fileSuffix+".tsv"
		else:
			outFileName = sys.argv[2]
	
		if narg < 4:
			username = input("Enter BRENDA username:")
			password = input("Enter BRENDA password:")
		else:
			username = sys.argv[3]
			password = sys.argv[4]
	
		password = hashlib.sha256(password.encode("utf-8")).hexdigest()
	
		main(organismName, outFileName, username, password)
	
