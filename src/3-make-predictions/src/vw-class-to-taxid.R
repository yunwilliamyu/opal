#Copyright (c) by respective owners including bioMerieux, and
#individual contributors Kevin Vervier and Pierre Mahe. All #rights reserved.  Released under a BSD (revised)
#license as described in the file LICENSE.

# read parameters from command line #
#-----------------------------------#
args = commandArgs(trailingOnly = TRUE)
input  = args[1]
dico = args[2]
output = args[3]


# read dictionary #
#-----------------#
dico = read.table(dico)

# read predictions #
#------------------#
preds = read.table(input)$V1

# convert #
#---------#
ind.match = match(preds, dico$V2)
if(sum(is.na(ind.match)) > 0){
	stop("ERROR : class(es) not found in dictionary")
}
preds.taxid = dico$V1[ind.match]

# save #
#------#
write.table(preds.taxid, file = output, quote = F, row.names = F, col.names = F)



