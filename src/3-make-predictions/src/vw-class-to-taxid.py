import sys

inputdir = sys.argv[1]
dicodir = sys.argv[2]
outputdir = sys.argv[3]

dico = {}
with open(dicodir, "r") as fin:
	for line in fin:
		txid, vwid = line.strip().split()[:2]
		dico[vwid] = txid
predout = open(outputdir, "w")
with open(inputdir, "r") as fin:
	for line in fin:
		predout.write("%s\n"%(dico[line.strip()]))
predout.close()