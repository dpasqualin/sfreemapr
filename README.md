Package Requirements for Production
====================
sudo apt-get install r-base-core

Package Requirements for Development
====================
sudo apt-get install r-base-core texlive-full

R Dependencies
==============
# Run within R
install.packages("phytools")

Build and Install
=================
R CMD build sfreemap && R CMD INSTALL sfreemap

Example
=======

require(devtools)
require(phytools)
require(sfreemap)
tree<-pbtree(n=100,scale=1)
Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-letters[1:nrow(Q)]
tree<-sim.history(tree,Q,anc="a")

sm<-sfreemap.map(tree, tree$states, Q='empirical')
sfreemap.pie_plot(sm)
