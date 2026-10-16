
library("rgl")
library("adegenet")
library("adegraphics")
library("treespace")
library("ggplot2")
library("ape")
library("tree")
library("phangorn")
library(parallel)
library(proxy)

get_dist_matrix = function(single_gene) {
  fasta_file = paste("~/Output_Nov5ccir_2rpt/pt06_Nov5_431sample_panaroo_onedir_core.99/aligned_gene_sequences/", single_gene, ".aln.fas", sep="")
  sequences <- read.dna(fasta_file, format = "fasta")
  print(sequences)
  dist_matrix <- dist.dna(sequences, model = "TN93")
  return(dist_matrix)
}

get_nj_tree = function(single_gene, tip_label_to_match=NULL) {
  fasta_file = paste("~/Output_Nov5ccir_2rpt/pt06_Nov5_431sample_panaroo_onedir_core.99/aligned_gene_sequences/", single_gene, ".aln.fas", sep="")
  sequences <- read.dna(fasta_file, format = "fasta")
  # print(sequences)
  dist_matrix <- dist.dna(sequences, model = "TN93")
  print(dim(as.matrix(dist_matrix)))
  # print(dist_matrix)
  length(dist_matrix)
  if( any(is.na( dist_matrix)) ) {
    return(NULL)
  }
  # Generate the Neighbor-Joining tree
  nj_tree <- njs(dist_matrix)
  nj_tree$tip.label <- unlist(lapply(strsplit(nj_tree$tip.label, split=";"), FUN=function(x)x[1]))
  nj_tree$tip.label <- sub("_R_", "", nj_tree$tip.label)
  cat("length(nj_tree$tip.label)=", length(nj_tree$tip.label), "\n")
  if( ! is.null(tip_label_to_match) ) {
    if(length(nj_tree$tip.label) != length(tip_label_to_match)) {
      return(NULL)
    }
    if( ! setequal(nj_tree$tip.label, tip_label_to_match) ) {
      return(NULL)
    }
  }
  # nj_tree$tip.label = c("")
  # plot(nj_tree, main = "NJ Tree for Single Gene")
  print(is.rooted(nj_tree))
  # midpoint.root(nj_tree)
  # root(nj_tree, method = "midpoint")
  # if(!is.null(outgroup)) {
  #   # nj_tree = root(nj_tree, outgroup = outgroup, resolve.root = FALSE, method = "default")
  #   nj_tree = root(nj_tree, outgroup = outgroup, method = "midpoint")
  # }
  nj_tree = midpoint(nj_tree)
  print(is.rooted(nj_tree))
  return(nj_tree)
}

# root(nj_tree, outgroup=c(100), method = "midpoint")
# midpoint(nj_tree)
# res = root_at_midpoint(nj_tree,
#                   split_edge      = TRUE,
#                   update_indices  = TRUE,
#                   as_edge_counts  = FALSE,
#                   is_rooted       = FALSE)
# root_tree = get_nj_tree(gene_names[10], outgroup=NULL)
# new_tree = get_nj_tree(gene_names[1], outgroup=root_tree)

do_contour = function(tree_list, cnt) {
  # my_mtree = tree_list
  # my_tspace <- treespace(my_mtree, method="treeVec", nf=3)
  my_tspace <- treespace(tree_list, method="treeVec", nf=3)
  class(my_tspace)
  pdf(paste("_", cnt, ".contour.pdf", sep=""))
  plot(my_tspace$pco$li[,1], my_tspace$pco$li[,2], main=paste("cnt=", cnt))
  my_tspace_df <- my_tspace$pco$li
  my_tspace_plot <- ggplot(my_tspace_df, aes(x=A1, y=A2)) # create plot
  my_tspace_plot2 <- my_tspace_plot + geom_density2d(colour="gray80") + # contour lines
  geom_point(size=6, shape=1, colour="gray50") + # grey edges
  geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
  xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background
  print(my_tspace_plot2)
  dev.off()
}

#NOT_WORK # Read in the phylogenetic tree (Newick format)
#NOT_WORK tree <- read.tree("aligned_genes.treefile")
#NOT_WORK
#NOT_WORK # Optionally, you can read in the alignment matrix (if you want to analyze the sequence data directly)
#NOT_WORK alignment <- read.dna("aligned_genes.fasta", format = "fasta")
#NOT_WORK
#NOT_WORK # Perform Treespace analysis on the phylogenetic tree
#NOT_WORK ts <- treespace(tree)
#NOT_WORK
#NOT_WORK # Plot the results
#NOT_WORK plot(ts)
#NOT_WORK #---------------------------------------------------------------------------------
#NOT_WORK
#NOT_WORK # a different approach
#NOT_WORK # Compute the distance matrix from the alignment
#NOT_WORK # dist_matrix <- dist.dna(alignment, model = "K80")  # Using Kimura 2-parameter model, for example
#NOT_WORK dist_matrix = get_dist_matrix(gene_names[1])
#NOT_WORK
#NOT_WORK # Perform Treespace analysis on the distance matrix
#NOT_WORK ts <- treespace(dist_matrix)
#NOT_WORK
#NOT_WORK # Plot the results
#NOT_WORK plot(ts)

# gene_names = scan("~/Output_Nov5ccir_2rpt/pt06_Nov5_431sample_panaroo_onedir_core.99/_cnt37.gene_names.txt", what="")
# gene_names = scan("~/Output_Nov5ccir_2rpt/pt06_Nov5_431sample_panaroo_onedir_core.99/_all_core_gene.txt", what="")
# gene_names = scan("~/Output_Nov5ccir_2rpt/pt06_Nov5_431sample_panaroo_onedir_core.99/_cnt1745.gene_names.txt", what="")
gene_names = scan("~/Output_Nov5ccir_2rpt/pt06_Nov5_431sample_panaroo_onedir_core.99/_largetst_cnt.gene_names.txt", what="")
# gene_name = gene_names[92] gene_names = gene_names[-92] gene_names = gene_names[-106]
# gene_names = gene_names[1:130]
tree_list = list()
tree_names = c()
cnt = 0
tip_label_to_match = NULL
for(i in seq(length(gene_names))) {
  cat("i=", i, "\n")
  gene_name = gene_names[i]
  res = get_nj_tree(gene_name, tip_label_to_match)
  # ape::write.tree(res, file=paste("_", i, "_nj.treefile", sep=""))
  if(! is.null(res)) {
     cnt = cnt + 1
     tree_list[[cnt]] = res
     tree_names[cnt] = gene_name
     if( is.null(tip_label_to_match) ) {
       tip_label_to_match = res$tip.label
     }
  }
  cat("cnt=", cnt, "\n")
  # if( cnt %% 150 == 0 ) {
  #   save(tree_list, file = paste("_", cnt, "_tree_list.RData", sep=""))
  #   # do_contour(tree_list, cnt)
  # }
}
class(tree_list) = "multiPhylo"
print(class(tree_list))
names(tree_list) <- tree_names
save(tree_list, file = paste("", cnt, "_tree_list.RData", sep=""))
do_contour(tree_list, cnt)
length(tree_list)

# tree1 = get_nj_tree(single_gene="zwf")
# tree2 = get_nj_tree(single_gene="yxeI")
# tree3 = get_nj_tree(single_gene="zapA")
# tree_list[[1]]$tip.label
# tree_list[[2]]$tip.label
# tree_list[[3]]$tip.label

# Error in storage.mode(x1) <- "double" :
#   'list' object cannot be coerced to type 'double'
# Calls: do_contour -> treespace -> as.dist -> rdist
# Execution halted

print("End: single_gene_tree.R")

# Note: see single_gene_tree.sh

