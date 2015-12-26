# THIS PACKAGE IS FOR TESTING ONLY

This is a first version of the *sfreemap* package, written entirely in R,
and used only for comparison between versions in my thesis.


### Package Requirements for Production

You need to have R installed on your system. If you are using a debian/ubuntu based distribution, just type the following command in a terminal.

`sudo apt-get install r-base-core`

### Package Requirements for Development

`sudo apt-get install r-base-core texlive-full`

###  R Dependencies

### Install

Before installing *sfreemapr* make sure you have the package *phytools*:

```
install.packages('phytools')
install.packages('devtools')
install_github('dpasqualin/sfreemapr')
```

If you have troubles installing the `devtools` package, try downloading
`sfreemapr` and then building and installing it using the following commands:

```
git clone https://github.com/dpasqualin/sfreemapr.git
R CMD check sfreemapr && R CMD build sfreemapr && R CMD INSTALL sfreemapr
```

If you choose to install using the command above, the documentation will be
available in the directory `sfreemapr.Rcheck`.


### Example

```
require(sfreemapr) # load package
tree <- pbtree(n=100,scale=1) # create a tree with 100 taxa
Q <- matrix(c(-1,1,1,-1),2,2) # create a transition rate matrix
rownames(Q)<-colnames(Q)<-letters[1:nrow(Q)] # give name to the states
tree <- sim.history(tree,Q,anc="a") # simulate a history for the tree

# estimate the history
sm <- sfreemapr.map(tree, tree$states, Q='empirical')
sfreemapr.pie_plot(sm) # plot the result
sfreemapr.describe(sm) # numerical summary of the result
```
