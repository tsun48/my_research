library(ape)
library(treespace)
load("1361_tree_list.RData")
print(length(tree_list))
print(names(tree_list))
class(tree_list) = "multiPhylo"
ts_res <- treespace(tree_list, nf=2) # default method= treeVec for Kendall Colijn metric vector (for rooted tree)
save(ts_res, file = paste("ts_res.RData", sep=""))

# cd ~/Output_Nov5ccir_2rpt/pt09_treespace_contour_plot 
# R --no-save < ~/Summer2024/ccirwork/save_ts_res.R
system("rsync -avz ts_res.RData remote_host:Output_Nov5ccir_2rpt/pt09_treespace_contour_plot/")


print("End: save_ts_res.R")


