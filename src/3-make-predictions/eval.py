import pandas as pd
import numpy as np
import argparse
from sklearn.metrics import precision_score, recall_score

parser = argparse.ArgumentParser()
parser.add_argument('-p', nargs=1)
parser.add_argument('-r', nargs=1)
args = parser.parse_args()
predfile = args.p[0]
reffile = args.r[0]

root = "../.."
#predfile = root + "/output/" + DB + "/" + ID + "/3-make-predictions/test.fragments."+DB+"-db.preds.taxid"
# tax id of generated test fragments
#reffile = root + "/output/" + DB + "/" + ID + "/1-generate-test-datasets/test.fragments.taxid"
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
       
