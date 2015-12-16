# Input
#   tree    a phylogenetic tree as an object of class "phylo" (from package
#           ape)
sfreemap.map <- function(tree, tip_states, Q, ...) {

    # Should this program run in parallel?
    parallel <- TRUE
    if (hasArg(parallel)) {
        parallel <- list(...)$parallel
    }

    # tree sanity check
    if (inherits(tree, "multiPhylo")) {
        # For Just call the same program multiple times...
        if (parallel == TRUE) {
            cores <- detectCores()
            mtrees <- mclapply(tree, sfreemap.map, tip_states, Q, ..., mc.cores=cores)
        } else {
            mtrees <- lapply(tree, sfreemap.map, tip_states, Q, ...)
        }

        # When Q=mcmc we will have length(trees)*n_simulation trees at the end
        # of the execution. Instead of lists of multPhylo objects we want to
        # return one multiPhylo object with all trees.
        if (Q == 'mcmc') {
            mtrees <- c(mapply(c, mtrees))
        }

        class(mtrees) <- "multiPhylo"
        return(mtrees)

    } else if (!inherits(tree, "phylo")) {
        stop("'tree' should be an object of class 'phylo'")
    } else if (!is.rooted(tree)) {
        stop("'tree' must be rooted")
    }

    # Defining the prior distribution for the root node of the tree,
    # also known as "pi"
    prior <- "equal"
    if (hasArg(prior)) {
        prior <- list(...)$prior
    }

    # tol gives the tolerance for zero elements in Q.
    # (Elements less then tol will be reset to tol)
    tol <- 1e-8
    if (hasArg(tol)) {
        tol <- list(...)$tol
    }

    # prior a list containing alpha and beta parameters for the gamma
    # prior distribution on the transition rates in Q. Note that alpha
    # and beta can be single values or vectors, if different prior
    # are desired for each value in Q
    gamma_prior <- list(alpha=1, beta=1, use.empirical=FALSE)
    if (hasArg(gamma_prior)) {
        pr <- list(...)$gamma_prior
        gamma_prior[names(pr)] <- pr
    }

    # burn_in for the MCMC
    burn_in <- 1000
    if (hasArg(burn_in)) {
        burn_in <- list(...)$burn_in
    }

    # sample_freq for the MCMC
    sample_freq <- 100
    if (hasArg(sample_freq)) {
        sample_freq <- list(...)$sample_freq
    }

    # number os simulations for the MCMC
    n_simulations <- 100
    if (hasArg(n_simulations)) {
        n_simulations <- list(...)$n_simulations
    }

    # For now only "symmetrical" model is accepted
    model <- "SYM"
    #if (hasArg(model)) {
    #    model <- list(...)$model
    #}

    # A single numeric value or a vector containing the (normal)
    # sampling variances for the MCMC
    vQ <- 0.1
    if (hasArg(vQ)) {
        vQ <- list(...)$vQ
    }

    # Define the tip states as a matrix
    if (!is.matrix(tip_states)) {
        tip_states <- build_states_matrix(tree, tip_states)
    }

    # Defining Q
    if (is.character(Q) && (Q == "empirical")) {
        # Phytools would replicate this result nsim times, but for now
        # we will return just one result.
        QP <- Q_empirical(tree, tip_states, prior, model, tol)
    } else if (is.character(Q) && (Q == "mcmc")) {
        # TODO: This function will generate many Qs. We have to decide how to
        # deal with it. Maybe just run the program for every Q?
        QP <- Q_mcmc(tree, tip_states, model, prior, gamma_prior, tol, burn_in
                     , sample_freq, vQ, n_simulations)
        # Call sfreemap.map for each {Q,prior} returned by the mcmc simulation
        params <- list(...)
        params$tree <- tree
        params$tip_states <- tip_states
        res <- function(QP) {
            params$Q <- QP$Q
            params$prior <- QP$prior
            return (do.call(sfreemap.map, params))
        }

        if (parallel == TRUE) {
            mtrees <- mclapply(QP, res, mc.cores=detectCores())
        } else {
            mtrees <- lapply(QP, res)
        }
        class(mtrees) <- "multiPhylo"
        return(mtrees)

    } else if (is.matrix(Q)) {
        # Phytools would replicate this result nsim times, but for now
        # we will return just one result.
        QP <- Q_matrix(tree, tip_states, Q, model, prior, tol)
    } else {
        stop("Unrecognized format for 'Q'")
    }

    # Set the final value
    Q <- QP$Q
    prior <- QP$prior
    logL <- QP$logL

    # Vector with rewards
    rewards <- rep(1,nrow(Q))
    if (hasArg(rewards)) {
        rewards <- list(...)$rewards
        if (length(rewards) != nrow(Q)) {
            stop("The rewards vector should represent the states")
        }
    }
    names(rewards) <- colnames(Q)

    # Defining the transitions of interest. By default, all transitions.
    QL <- matrix(1, nrow=nrow(Q), ncol=ncol(Q))
    diag(QL) <- 0
    if (hasArg(QL)) {
        QL <- list(...)$QL
        if (!all(dim(Q) == dim(QL))) {
            stop("QL must have same dimensions as Q")
        }
    }
    rownames(QL) <- colnames(QL) <- rownames(Q)

    # Acquire more info about the tree.
    tree_extra <- list(
        states = tip_states
        , n_states = nrow(Q)
        , n_edges = length(tree$edge.length)
        , n_tips = nrow(tip_states)
        , n_nodes = nrow(tip_states) + tree$Nnode
        , rewards = rewards
    )

    # Reorder the tree so the root is the first row of the matrix.
    # We save the original order to make sure we have the result
    # in same order of the tree;
    tree <- reorder(tree, 'pruningwise')

    # Step 1
    # Compute Eigen values and Eigen vectors for the transition rate matrix
    # Q, assuming Q is symmetrical
    Q_eigen <- eigen(Q, TRUE, only.values = FALSE)
    # The inverse of Q_eigen vectors
    Q_eigen[['vectors_inv']] <- solve(Q_eigen$vectors)

    # Step 2
    # Compute P(tp), the transistion probability, for each edge length t of T
    # Q = U X diag(d1,...,dm) X U**-1
    # U are the eigenvectors of Q
    # d1,...,dm are the eigenvalues of Q
    # diag(d1,...,dm) is a diagonal matrix with d1,...,dm on it's main
    # diagonal
    # For now I'm doing this just where it is needed, in the
    # fractional_likelihood calculation. I don't know yet if I'm going
    # to need it somewhere else. It might be useful to compute it here
    # and use it in all different places if that's the case.
    MAP <- list()

    MAP[['Q']] <- Q
    MAP[['prior']] <- prior

    # Transistion probabilities
    MAP[['tp']] <- transition_probabilities(Q_eigen, tree$edge.length)

    # Step 4 and 5
    # Traverse the tree once and calculate Fu and Sb for each node u and
    # each edge b;
    # Compute the data likelihood Pr(D) as the dot product of Froot and root
    # distribution pi.
    MAP[['fl']] <- fractional_likelihoods(tree, tree_extra, Q, Q_eigen
                                          , prior, MAP$tp, tol)


    # FIXME: for some unknown reason if I remove this line I get a 'out of
    # bounds' error in posterior_restricted_moment(). This line doesn't need to
    # exist
    MAP[['h']] <- list()

    # Build dwelling times
    tree[['mapped.edge']] <- tree[['mapped.edge.lmt']] <- tname <- NULL
    states <- rownames(QL)
    n_states <- tree_extra$n_states
    multiplier <- matrix(0, nrow=n_states, ncol=n_states)

    for (i in 1:n_states) {
        for (j in 1:n_states) {

            if (i == j) {
                value <- rewards[i] # dwelling times
            } else {
                value <- QL[i,j] * Q[i,j] # number of transitions
            }

            if (value == 0) {
                next
            }
            multiplier[i,j] <- value

            # Step 3
            h <- func_H(multiplier, Q_eigen, tree, tree_extra, tol)
            prm <- posterior_restricted_moment(h, tree, tree_extra, MAP)
            ev <- expected_value(tree, Q, MAP, prm)

            if (i == j) {
                # dwelling times
                tree[['mapped.edge']] <- cbind(tree[['mapped.edge']], ev)
            } else {
                # number of transitions
                state_from <- states[i]
                state_to <- states[j]
                tname <- c(tname, paste(state_from, state_to, sep=','))
                tree[['mapped.edge.lmt']] <- cbind(tree[['mapped.edge.lmt']], ev)
            }

            multiplier[i,j] <- 0
        }
    }

    colnames(tree[['mapped.edge']]) <- names(rewards)
    colnames(tree[['mapped.edge.lmt']]) <- tname

    # Return the tree in the original order
    return (sfreemap.reorder(tree, 'cladewise'))
}

# The final answer!
expected_value <- function(tree, Q, map, prm) {

    likelihood <- map[['fl']][['L']]

    ret <- prm / likelihood

    # the rownames of the mapped objects
    names(ret) <- paste(tree$edge[,1], ",", tree$edge[,2], sep="")

    return (ret)
}

# Compute the expected value for labelled markov transitions and
# evolutionary reward, given the fractional and directional likelihoods
posterior_restricted_moment <- function(multiplier, tree, tree_extra, map) {
    # Allocate vector of size tree_extra$n_edges. Index i will hold the expected value
    # of the ith edge of the tree
    ret <- rep(0, tree_extra$n_edges)

    fl <- map[['fl']] # fractional likelihoods
    F <- fl[['F']]
    G <- fl[['G']]
    S <- fl[['S']]

    # These tree lines below are responsible to generate a matrix similar to
    # tree$edge, but with an additional columns representing the 'brother'
    # of the node.
    get_brother <- function(x) return(tree$edge[c(x+1,x),2])
    in_pairs <- seq(1, nrow(tree$edge), 2)
    edge <- cbind(tree$edge, sapply(sapply(in_pairs, get_brother), rbind))

    for (i in 1:tree_extra$n_edges) {
        p <- edge[i,1] # parent node
        c <- edge[i,2] # child node
        b <- edge[i,3] # brother node

        partial <- multiplier[,,i]

        # Equation 3.8 in the paper
        for (j in 1:tree_extra$n_states) {
            gs <- G[p,j] * S[b,j] * F[c,]
            ret[i] <- ret[i] + sum(gs * partial[j,])
        }
    }

    return (ret)
}

# The vector of forward, often called partial or fractional likelihood.
fractional_likelihoods <- function(tree, tree_extra, Q, Q_eigen, prior, Tp, tol) {

    # F is the vector of forward, often called partial or fractional,
    # likelihoods at node u. Element F ui is the probability of the
    # observed data at only the tips that descend from the node u,
    # given that the state of u is i.
    F <- matrix(0, tree_extra$n_nodes, tree_extra$n_states)
    # Directional likelihoods
    S <- matrix(0, tree_extra$n_nodes, tree_extra$n_states)
    # G is the probability of observing state i at node
    # u together with other tip states on the subtree of t
    # obtained by removing all lineages downstream of node u.
    G <- matrix(0, tree_extra$n_nodes, tree_extra$n_states)

    # Set col names, just to make it easy to debug
    states <- colnames(F) <- colnames(S) <- colnames(G) <- colnames(Q)

    # Init F matrix with tip values
    # As stated in the article, when value is ambiguous, we set all ambiguous
    # values to 1 (100% chance)
    F[1:tree_extra$n_tips,] <- tree_extra$states

    # Compute F
    for (e in seq(1, tree_extra$n_edges, 2)) {
        u <- tree$edge[e,1]        # parent node
        right <- tree$edge[e,2]    # right node
        left <- tree$edge[e+1,2]   # left node

        # NOTE: this loop can possibly be vectorized...
        tright <- Tp[,,e]     # transistion probability for edge p -> right
        tleft <- Tp[,,e+1]    # transistion probability for edge p -> left
        for (i in 1:tree_extra$n_states) {
            S[right,i] <- sum(F[right,] * tright[i,])
            S[left,i] <- sum(F[left,] * tleft[i,])
        }

        # When values are smaller than tol set them to tol, otherwise the
        # likelihood will be zero, we will have division by zero and the world
        # will end
        F[u,] <- S[right,] * S[left,]
    }

    # Get the root node number
    # Remember that the tree is ordered pruninwise
    root <- tail(tree$edge, n=1)[1,1]

    # Compute the likelihood of the tree
    L <- sum(F[root,]*prior)

    # Compute G using F and S
    # Article: "Computational Advances in Maximum Likelihood Methods for
    # Molecular Phylogeny", authors Eric E. Schadt, Janet S. Sinsheimer and
    # Kenneth Lange
    G[root,] <- prior
    for (e in seq(tree_extra$n_edges, 1, -2)) {
        p <- tree$edge[e,1]       # parent node
        left <- tree$edge[e,2]    # left node
        right <- tree$edge[e-1,2] # right node

        # TODO: Vectorize this...
        tleft <- Tp[,,e]
        tright <- Tp[,,e-1]
        for (i in 1:tree_extra$n_states) {
            G[left,i] <- sum(G[p,]*S[right,]*tleft[i,])
            G[right,i] <- sum(G[p,]*S[left,]*tright[i,])
        }
    }

    return (list(F=F, G=G, S=S, L=L))
}

# Input
# Q_eigen   the eigen deposition of Q
# t         a edge length
#
# output    A tridimentional matrix with the transition probabilities P
#           for every the edge length t of the tree.
transition_probabilities <- function(Q_eigen, edges) {

    # P(t) = U X diag(e**d1t, ... e**dmt) X U**-1
    # t is an arbitrary edge length
    # exp(n) is nth power of e (euler number)
    n_edges <- length(edges)
    P <- array(0, dim = c(dim(Q_eigen$vectors), n_edges))

    # Initialize d
    d <- diag(Q_eigen$values)

    # Iterate over branch length
    # NOTE: vectorize this...
    for (i in 1:n_edges) {
        diag(d) <- exp(Q_eigen$values * edges[i])
        P[,,i] <- Q_eigen$vectors %*% d %*% Q_eigen$vectors_inv
    }

    return(P)
}

# The expected number of labelled markov transitions and expected markov rewards
# TODO: find a better name for this function =P
func_H <- function(multiplier, Q_eigen, tree, tree_extra, tol) {

    ret <- array(0, dim = c(tree_extra$n_states, tree_extra$n_states, tree_extra$n_edges))

    # Just an alias, to make it easier to read the formulas below
    d <- Q_eigen$values

    # 1. Get the transitions that happened in a edge
    # TODO: This matrix has only one entry different from zero, there
    # might be a faster way to calculate Si and Sj without multiplying
    # the entire matrix
    gen_Si <- function(i) {
        E <- matrix(0, tree_extra$n_states, tree_extra$n_states)
        E[i,i] <- 1
        return (Q_eigen$vectors %*% E %*% Q_eigen$vectors_inv)
    }

    # This is part of the second inner loop. We broght it outside
    # because... it's R, and you have to do weird stuff to improve performance.
    Sij_all <- lapply(1:tree_extra$n_states, gen_Si)
    S_partial <- lapply(1:tree_extra$n_states, function(i) Sij_all[[i]] %*% multiplier)

    # A matrix where each line corresponds to exp(t*d) for the nth edge
    # This is gonna be used inside build_Iij function. It was there before,
    # But doing this way makes it faster. Just R stuff...
    Ib_all <- do.call(rbind, lapply(tree$edge.length, function(t) exp(d*t)))

    for (b in 1:tree_extra$n_edges) {
        t <- tree$edge.length[b]

        # NOTE: Maybe initialize outside the loop and just zero it inside?
        partial <- matrix(0, tree_extra$n_states, tree_extra$n_states)

        Ib <- Ib_all[b,]

        for (i in 1:tree_extra$n_states) {
            # 3. S, which uses the Q_eigen and E, a matrix with zero entries
            # except for the iith element, which is 1
            # Set new Ei
            Si_partial = S_partial[[i]]

            for (j in 1:tree_extra$n_states) {
                # Set new Ej
                Sj <- Sij_all[[j]]

                # 4. Iij, which uses the edge length t and the eigen values di,dj
                # TODO: This can be moved to outside both for loops
                Iij <- build_Iij(Ib, d, i, j, t, tol)

                partial <- partial + ((Si_partial %*% Sj) * Iij)
            }
        }

        ret[,,b] <- partial
    }

    return (ret)
}

# Maybe its better to build the entire matrix with all
# possible values once and just use its values afterwards
build_Iij <- function(Ib, d, i, j, t, tol) {

    # NOTE: We are comparing two floats here,
    # tolerance got from (Paula Tataru, 2011)
    # make.simmap (phytools) accepts this tolerance as a
    # parameter, should we have it as a parameter too?
    # NOTE: NEVER use all.equal() when not really necessary, it's absurdly
    # slow..
    diff <- d[i] - d[j]
    if (abs(diff) <= tol) {
        return ( t * Ib[i] )
    } else {
        return ( (Ib[i] - Ib[j]) / diff )
    }
}
