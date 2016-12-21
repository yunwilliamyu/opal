# specify model parameters 
DB = "A1"
NBATCHES = "2"

computePerf = function(ref,pred){
	pred.correct = ref == pred
	acc.micro = mean(pred.correct)
	tmp = split(pred.correct, ref)
        acc.micro = round( 100*mean(pred.correct), digits = 2)
	acc.species = round( 100*sapply(tmp, mean), digits = 2)
	acc.macro = round(mean(acc.species), digits = 2)
	acc.median = round(median(acc.species), digits = 2)
	return(list("acc.micro"=acc.micro,"acc.macro"=acc.macro,"acc.median"=acc.median,"acc.species"=acc.species))
}


perf = list()

# read results fragments 
# read reference labels
ref = read.table("../../1-generate-test-datasets/output/test.fragments.taxid")$V1
# read predictions
pred = read.table(paste("../output/test.fragments.", DB, "-db.preds.taxid", sep = ""))$V1
# compute perf
perf.tmp = computePerf(ref, pred)
perf[["fragments"]] = perf.tmp
perf.tmp["acc.micro"]
perf.tmp["acc.macro"]
perf.tmp["acc.median"]
