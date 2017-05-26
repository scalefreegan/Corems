get.mast <- function (gene, window = 125, shift = 75, e.value.cutoff = Inf,
    p.value.cutoff = 1e-06, op.shift = F, include.bad = F, verbose = T,
    biclust.filter = NULL, motif.filter = NULL, count.all = F) {
    if (is.character(gene)) {
        coo <- e$get.gene.coords(gene, op.shift = op.shift)
        if (is.null(coo))
            return(NULL)
        print(coo)
        st.st <- c(coo$start_pos, coo$end_pos)
        if (coo$strand == "R")
            st.st <- rev(st.st)
        st.st <- st.st[1] + c(-window, +window) + (if (coo$strand ==
            "R")
            +shift
        else -shift)
        chr <- as.character(coo$contig)
        names(st.st)[1] <- chr
    } else{
      print("Input must be gene")
      return(NULL)
    }
    if (!exists("pssm.scans"))
        pssm.scans <- get.pssm.scans()
    if (p.value.cutoff < max(pssm.scans$pvals, na.rm = T))
        pssm.scans <- pssm.scans[pvals <= p.value.cutoff]
    if (!is.data.table(pssm.scans)) {
        pssm.scans <- as.data.table(pssm.scans)
        gc()
        setkey(pssm.scans, bic, mots, gene, posns)
    }
    if (st.st[1] < 1)
        st.st[1] <- 1
    chr <- names(st.st)[1]
    if (is.null(chr)) {
        chr <- names(which(sapply(e$genome.info$genome.seqs,
            nchar) == max(sapply(e$genome.info$genome.seqs, nchar))))
        names(st.st)[1] <- chr
    }
    scans <- pssm.scans[gene == chr & posns %betw% (st.st +
        c(-100, 100))]
    scans <- unique(scans)
    motifs <- unique(paste("MOT", scans$bic, abs(scans$mots),
        sep = "_"))
    if (length(motifs) <= 0)
        stop("No motifs pass criteria (1)!")
    if (verbose)
        cat(length(motifs), "motifs.\n")
    if (!is.null(biclust.filter)) {
        motifs <- motifs[motifs %chin% unlist(get.motifs(biclust = biclust.filter))]
        if (length(motifs) <= 0)
            stop("No motifs pass criteria! (2)")
    }
    if (!is.null(motif.filter)) {
        motifs <- motifs[motifs %chin% motif.filter]
        if (length(motifs) <= 0)
            stop("No motifs pass criteria! (3)")
    }
    if (!include.bad && exists("bad.clusts")) {
        bad.ms <- unique(unlist(get.motifs(motif.clust = bad.clusts,
            expand = F)))
        motifs <- motifs[!motifs %chin% bad.ms]
        if (length(motifs) <= 0)
            stop("No motifs pass criteria! (4)")
    }
    if (!is.infinite(e.value.cutoff) && !is.na(e.value.cutoff)) {
        minfo <- get.motif.info(motifs = motifs)
        e.vals <- do.call(c, lapply(minfo, function(tmp) {
            if (is.null(tmp))
                return(NA)
            return(tmp$e.value)
        }))
        motifs <- motifs[e.vals <= e.value.cutoff]
        if (length(motifs) <= 0)
            stop("No motifs pass criteria! (5)")
    }
    if (!include.bad) {
        if (exists("coding.fracs")) {
            frac.in.coding <- coding.fracs$all.fracs[motifs]
        }
        else {
            coding.fracs <- get.motif.coding.fracs(motifs, verbose = T)
            frac.in.coding <- coding.fracs$all.fracs
        }
        motifs.orig <- motifs
        motifs <- motifs[!is.na(frac.in.coding) & frac.in.coding <
            coding.fracs$mean.fracs - 0.01]
        if (length(motifs) <= 0)
            stop("No motifs pass criteria! (6)")
        rm(coding.seqs, scans, in.coding)
    }
    if (verbose)
        cat(length(motifs), "motifs remain.\n")
    mots <- strsplit(gsub("MOT_", "", motifs), "_")
    bi <- as.integer(sapply(mots, "[", 1))
    mo <- as.integer(sapply(mots, "[", 2))
    scans <- pssm.scans[J(c(bi, bi), c(mo, -mo), chr), allow.cart = T]
    scans <- scans[!is.na(scans$posns), ]
    scans <- scans[scans$posns %betw% (st.st + c(-500, 500)),
        ]
    setkey(scans, "bic", "mots")
    getEntropy <- function(pssm) {
        pssm[pssm == 0] <- 1e-05
        entropy <- apply(pssm, 1, function(i) -sum(i * log2(i)))
        return(entropy)
    }
    seq <- substr(e$genome.info$genome.seqs[names(st.st)[1]],
        st.st[1], st.st[2])
    mat <- matrix(0, nrow = diff(st.st) + 1, ncol = 4)
    rownames(mat) <- as.character(st.st[1]:st.st[2])
    colnames(mat) <- e$col.let
    scans <- scans[scans$posns %betw% (st.st + c(-100, 100))]
    if (exists("motif.clusts")) {
        mcs <- get.motif.clusters(motif = motifs)
        mot.tab <- sort(table(unlist(mcs)))
        mot.tab <- rev(mot.tab[mot.tab > 2])
        print(mot.tab)
    }
    else {
        mot.tab <- character()
    }
    if (length(mot.tab) > 10)
        mot.tab <- mot.tab[1:9]
    if (count.all)
        mot.tab <- c(ALL = 0, mot.tab)
    counts <- matrix(0, nrow = nrow(mat), ncol = length(mot.tab))
    colnames(counts) <- names(mot.tab)
    rownames(counts) <- rownames(mat)
    if (nrow(scans) > 0) {
        bics <- unique(scans$bic)
        for (k in bics) {
            if (verbose) {
                wh <- which(bics == k)
                if (wh%%100 == 1)
                  cat(k, wh, length(bics), "\n")
            }
            sc <- scans[bic == k]
            mots <- unique(abs(sc$mots))
            for (m in mots) {
                width <- motif.widths[k, m]
                if (width <= 0)
                  next
                if (exists("mcs")) {
                  mc <- mcs[[paste("MOT", k, m, sep = "_")]]
                  mc <- mc[mc %chin% colnames(counts)]
                }
                pssm.orig <- get.motif.info(paste("MOT", k, m,
                  sep = "_"))[[1]]$pssm
                pssm.rev <- pssm.orig[nrow(pssm.orig):1, 4:1]
                entr <- getEntropy(pssm.orig)
                scale.e.orig <- (2 - entr)/2
                scale.e.rev <- rev(scale.e.orig)
                sc2 <- sc[abs(sc$mots) == m]
                inds.pssm <- 1:nrow(pssm.orig)
                for (i in 1:nrow(sc2)) {
                  mot <- sc2$mots[i]
                  if (sign(mot) == 1) {
                    pssm <- pssm.orig
                  }
                  else if (sign(mot) == -1) {
                    pssm <- pssm.rev
                  }
                  posn <- sc2$posns[i]
                  inds <- (posn - 1):(posn - 2 + width)
                  inds.1 <- inds.pssm[inds %betw% st.st]
                  if (length(inds.1) <= 0)
                    next
                  inds <- inds[inds %betw% st.st]
                  inds <- inds - st.st[1] + 1
                  mat[inds, ] <- mat[inds, ] + pssm[inds.1, ]
                  if (count.all)
                    counts[inds, "ALL"] <- counts[inds, "ALL"] +
                      1
                  if (exists("mc") && length(mc) > 0 && !is.null(mc)) {
                    counts[inds, mc] <- counts[inds, mc] + 1
                  }
                }
            }
        }
    }
    # return(invisible(list(motifs = motifs, scans = scans,
    #     mat = mat, counts = counts, st.st = st.st, mot.tab = mot.tab)))
    return(scans)
}
