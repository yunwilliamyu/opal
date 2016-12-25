import pandas as pd
import numpy as np
import argparse
from sklearn.metrics import precision_score, recall_score

LOCAL = False
if LOCAL:
	DB = "A1"
else:
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', nargs=1)
	args = parser.parse_args()
	DB = args.d[0]

predfile = "../output/test.fragments."+DB+"-db.preds.taxid"
# tax id of generated test fragments
reffile = "../../1-generate-test-datasets/output_"+DB+"/test.fragments.taxid"
#reffile = "../../1-generate-test-datasets/output/test.fragments.taxid"
# tax id of existing test fragmetns
#reffile = "../../../data/"+DB+"/test/"+DB+".test.taxid"
with open(predfile, "r") as fin:
	pred = fin.read().splitlines()
with open(reffile, "r") as fin:
	ref = fin.read().splitlines()

pred = map(int, pred)
ref = map(int, ref)
correct = np.equal(pred, ref)

perf = pd.DataFrame({"pred":pred, "ref":ref, "correct":correct})
tmp = perf.groupby("ref")
species = tmp["correct"].agg(np.mean)
micro = np.mean(correct)
macro = np.mean(species)
median = np.median(species)
print "micro = %.2f\nmacro = %.2f\nmedian = %.2f"\
        %(micro*100, macro*100, median*100)
       

#print "prec(macro)=%.6f"%(precision_score(ref, pred, average='macro')*100)
#print "prec(micro)=%.6f"%(precision_score(ref, pred, average='micro')*100)
#print "recall(macro)=%.6f"%(recall_score(ref, pred, average='macro')*100)
#print "recall(micro)=%.6f"%(recall_score(ref, pred, average='micro')*100)

