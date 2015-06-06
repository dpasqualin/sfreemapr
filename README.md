### Package Requirements for Production

You need to have R installed on your system. If you are using a debian/ubuntu based distribution, just type the following command in a terminal.

`sudo apt-get install r-base-core`

### Package Requirements for Development

`sudo apt-get install r-base-core texlive-full`

###  R Dependencies

### Install

Before installing *sfreemap* make sure you have the package *phytools*:

```
install.packages('phytools')
install.packages('devtools')
install_github('dpasqualin/sfreemap')
```

If you have troubles installing the `devtools` package, try downloading
`sfreemap` and then building and installing it using the following commands:

```
git clone https://github.com/dpasqualin/sfreemap.git
R CMD check sfreemap && R CMD build sfreemap && R CMD INSTALL sfreemap
```

If you choose to install using the command above, the documentation will be
available in the directory `sfreemap.Rcheck`.


### Example

```
require(sfreemap) # load package
tree <- pbtree(n=100,scale=1) # create a tree with 100 taxa
Q <- matrix(c(-1,1,1,-1),2,2) # create a transition rate matrix
rownames(Q)<-colnames(Q)<-letters[1:nrow(Q)] # give name to the states
tree <- sim.history(tree,Q,anc="a") # simulate a history for the tree

# estimate the history
sm <- sfreemap.map(tree, tree$states, Q='empirical')
sfreemap.pie_plot(sm) # plot the result
sfreemap.describe(sm) # numerical summary of the result
```
