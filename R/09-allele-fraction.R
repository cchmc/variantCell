# ---------------------------------------------------------------------------
# Allele-fraction testing.
#
# The package already tests genotype (findSNPsByGroup, presence/absence) and
# read depth (findDESNPs, transcript abundance). Neither tests the quantity
# AD/DP itself - the fraction of reads carrying the alternative base. That
# fraction is the natural statistic for anything where the variant is in the
# RNA rather than the DNA: A-to-I editing rate, allele-specific expression,
# mtDNA heteroplasmy.
#
# Two entry points, because they answer different questions at very different
# statistical cost:
#
#   computeAlleleFractionIndex()  aggregates AD/DP over a site set into ONE
#                                 number per pseudobulk. One test, so no
#                                 multiple-testing burden. This is the Alu
#                                 Editing Index construction, and at the group
#                                 sizes available in this cohort it is the only
#                                 one of the two that can reach significance.
#
#   findDEAlleleFraction()        tests each site separately. Useful on a
#                                 restricted site list (recoding sites, a
#                                 candidate panel) or in a large cohort. On a
#                                 genome-wide editing site set at n < 10 it is
#                                 arithmetically incapable of clearing
#                                 Bonferroni, and says so rather than returning
#                                 a table of hopeful p-values.
#
# The unit of observation is the SAMPLE, never the cell. Cells within a library
# share a genome, a capture and a sequencing run; treating 8,000 of them as
# independent observations is pseudoreplication and inflates significance by
# orders of magnitude. Cells are aggregated to pseudobulk first, and the test
# runs across pseudobulks.
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
#' Smallest attainable two-sided p for a paired test over n pairs, all
#' differences in the same direction. Like .min_attainable_p() for the unpaired
#' case this is a permutation-count ceiling: it depends only on n, so no depth,
#' cell count or choice of statistic can beat it.
.min_attainable_p_paired <- function(n) {
  if (is.na(n) || n < 1) return(NA_real_)
  min(1, 2 / (2^n))
}

#' @keywords internal
#' @noRd
#' Resolve a user site specification to row indices of the SNP matrices.
#' Accepts integer indices, rs IDs, or "chr:pos" strings; NULL means all sites.
.resolve_sites <- function(snp_database, sites) {
  n <- nrow(snp_database$snp_info)
  if (is.null(sites)) return(seq_len(n))
  if (is.numeric(sites)) {
    bad <- sites < 1 | sites > n
    if (any(bad)) stop("site indices out of range", call. = FALSE)
    return(as.integer(sites))
  }
  info <- snp_database$snp_info
  key  <- paste0(info$CHROM, ":", info$POS)
  idx  <- match(sites, key)
  if (!is.null(info$rs_id)) {
    miss <- is.na(idx)
    if (any(miss)) idx[miss] <- match(sites[miss], info$rs_id)
  }
  if (all(is.na(idx))) {
    stop("none of the requested sites were found in the database", call. = FALSE)
  }
  if (any(is.na(idx))) {
    warning(sprintf("%d of %d requested sites not found and dropped",
                    sum(is.na(idx)), length(sites)), call. = FALSE)
  }
  as.integer(idx[!is.na(idx)])
}

#' @keywords internal
#' @noRd
#' Sum AD and DP over cells within each level of `groups`, returning
#' sites x groups matrices. Vectorised via an indicator matrix; a loop over
#' 250k sites would dominate the runtime.
.aggregate_pseudobulk <- function(AD, DP, groups) {
  keep <- !is.na(groups)
  if (!any(keep)) stop("no cells remain after grouping", call. = FALSE)
  g <- droplevels(factor(groups[keep]))
  ind <- Matrix::sparse.model.matrix(~ 0 + g)
  colnames(ind) <- levels(g)
  list(AD = as.matrix(AD[, keep, drop = FALSE] %*% ind),
       DP = as.matrix(DP[, keep, drop = FALSE] %*% ind),
       n_cells = as.integer(table(g)))
}

#' @keywords internal
#' @noRd
#' Method-of-moments overdispersion (rho) for a beta-binomial, pooled across
#' sites. Returns 0 when the data are not overdispersed relative to binomial.
.estimate_rho <- function(ad, dp) {
  ok <- dp > 0
  if (sum(ok) < 3) return(0)
  p <- sum(ad[ok]) / sum(dp[ok])
  if (p <= 0 || p >= 1) return(0)
  obs <- sum((ad[ok] - dp[ok] * p)^2 / (dp[ok] * p * (1 - p)))
  df  <- sum(ok) - 1
  if (df < 1) return(0)
  phi <- obs / df
  n_bar <- mean(dp[ok])
  rho <- (phi - 1) / (n_bar - 1)
  max(0, min(rho, 0.99))
}

#' @keywords internal
#' @noRd
#' Two-group beta-binomial likelihood ratio test at a single site.
.betabinom_lrt <- function(ad1, dp1, ad2, dp2, rho) {
  ll <- function(ad, dp, p) {
    if (p <= 0) p <- 1e-9
    if (p >= 1) p <- 1 - 1e-9
    if (rho <= 0) return(sum(stats::dbinom(ad, dp, p, log = TRUE)))
    s <- (1 - rho) / rho
    a <- p * s; b <- (1 - p) * s
    sum(lbeta(ad + a, dp - ad + b) - lbeta(a, b) + lchoose(dp, ad))
  }
  tot1 <- sum(dp1); tot2 <- sum(dp2)
  if (tot1 == 0 || tot2 == 0) return(NA_real_)
  p0 <- (sum(ad1) + sum(ad2)) / (tot1 + tot2)
  p1 <- sum(ad1) / tot1
  p2 <- sum(ad2) / tot2
  stat <- 2 * ((ll(ad1, dp1, p1) + ll(ad2, dp2, p2)) -
               (ll(ad1, dp1, p0) + ll(ad2, dp2, p0)))
  if (!is.finite(stat) || stat < 0) return(NA_real_)
  stats::pchisq(stat, df = 1, lower.tail = FALSE)
}


#' @title computeAlleleFractionIndex: Aggregate Alt-Allele Fraction per Pseudobulk
#' @name computeAlleleFractionIndex
#'
#' @description
#' Collapses AD/DP over a set of sites into a single alt-fraction index per
#' pseudobulk, where a pseudobulk is one level of `split_by` (normally sample)
#' crossed with one level of `group_by` and optionally `within`.
#'
#' This is the Alu Editing Index construction generalised: one number per
#' pseudobulk rather than one test per site. Because it produces a single
#' statistic it carries no genome-wide multiple-testing burden, which is what
#' makes it usable at the group sizes typical of a transplant cohort.
#'
#' @param group_by Character. Metadata column defining the arms being compared.
#'   Defaults to the current project identity.
#' @param split_by Character. Metadata column defining the unit of observation.
#'   Defaults to "sample_id". Do not set this to a per-cell column - the whole
#'   point is that cells are not independent.
#' @param within Character or NULL. Optional further stratification, e.g. a cell
#'   type column, producing one index per sample x arm x cell type.
#' @param sites Optional site restriction: row indices, rs IDs, or "chr:pos".
#'   NULL uses all sites in the database.
#' @param min_dp_site Numeric. Minimum pseudobulk depth for a site to contribute.
#' @param min_cells Numeric. Minimum cells for a pseudobulk to be reported.
#' @param common_sites Logical. If TRUE (default), restricts every pseudobulk
#'   within a `split_by` unit to sites covered in all of that unit's arms, so an
#'   index difference cannot be produced by the two arms using different sites.
#' @param weighted Logical. If TRUE (default) the index is the read-weighted
#'   sum(AD)/sum(DP). If FALSE it is the unweighted mean of per-site fractions,
#'   which is dominated by low-coverage sites.
#' @param verbose Logical. Print progress.
#'
#' @return A data frame with one row per pseudobulk: split, group, within,
#'   n_cells, n_sites, total_ad, total_dp, index. Attributes carry the
#'   parameters used and the paired-test ceiling for the observed n.
#'
#' @examples
#' \dontrun{
#' project$setProjectIdentity("donor_type")
#' idx <- project$computeAlleleFractionIndex(within = "cell_type")
#' }
#' @export
variantCell$set("public", "computeAlleleFractionIndex", function(group_by = NULL,
                                                                 split_by = "sample_id",
                                                                 within = NULL,
                                                                 sites = NULL,
                                                                 min_dp_site = 10,
                                                                 min_cells = 20,
                                                                 common_sites = TRUE,
                                                                 weighted = TRUE,
                                                                 verbose = TRUE) {
  if (is.null(self$snp_database)) stop("No SNP database. Run buildSNPDatabase() first.", call. = FALSE)
  meta <- self$snp_database$cell_metadata
  if (is.null(group_by)) group_by <- self$current_project_ident
  if (is.null(group_by)) stop("No identity set. Pass group_by or run setProjectIdentity().", call. = FALSE)
  for (col in c(group_by, split_by, within)) {
    if (!col %in% colnames(meta)) stop(sprintf("Column '%s' not found in cell metadata", col), call. = FALSE)
  }

  site_idx <- .resolve_sites(self$snp_database, sites)
  AD <- self$snp_database$ad_matrix[site_idx, , drop = FALSE]
  DP <- self$snp_database$dp_matrix[site_idx, , drop = FALSE]

  # One pseudobulk per split x group (x within). The key is kept as a single
  # factor so the aggregation is a single sparse matrix product.
  parts <- list(meta[[split_by]], meta[[group_by]])
  if (!is.null(within)) parts <- c(parts, list(meta[[within]]))
  key <- do.call(paste, c(parts, sep = "\r"))
  key[Reduce(`|`, lapply(parts, is.na))] <- NA

  if (verbose) cat(sprintf("Aggregating %d sites over %d pseudobulks...\n",
                           length(site_idx), length(unique(stats::na.omit(key)))))
  agg <- .aggregate_pseudobulk(AD, DP, key)

  fields <- do.call(rbind, strsplit(colnames(agg$DP), "\r", fixed = TRUE))
  out <- data.frame(
    split   = fields[, 1],
    group   = fields[, 2],
    stringsAsFactors = FALSE
  )
  if (!is.null(within)) out$within <- fields[, 3]
  out$n_cells <- agg$n_cells

  # Restrict to sites callable in every arm of the same split unit, so the two
  # arms are compared on identical positions.
  usable <- agg$DP >= min_dp_site
  if (common_sites) {
    block <- if (is.null(within)) out$split else paste(out$split, out$within, sep = "\r")
    for (b in unique(block)) {
      cols <- which(block == b)
      if (length(cols) > 1) {
        keep <- Reduce(`&`, lapply(cols, function(j) usable[, j]))
        for (j in cols) usable[, j] <- keep
      }
    }
  }

  idx <- vapply(seq_len(ncol(agg$DP)), function(j) {
    ok <- usable[, j]
    if (!any(ok)) return(NA_real_)
    if (weighted) sum(agg$AD[ok, j]) / sum(agg$DP[ok, j])
    else mean(agg$AD[ok, j] / agg$DP[ok, j])
  }, numeric(1))

  out$n_sites  <- colSums(usable)
  out$total_ad <- vapply(seq_len(ncol(agg$DP)), function(j) sum(agg$AD[usable[, j], j]), numeric(1))
  out$total_dp <- vapply(seq_len(ncol(agg$DP)), function(j) sum(agg$DP[usable[, j], j]), numeric(1))
  out$index    <- idx

  out <- out[out$n_cells >= min_cells & !is.na(out$index), , drop = FALSE]
  out <- out[order(out$group, out$split), , drop = FALSE]
  rownames(out) <- NULL

  n_pairs <- length(unique(out$split))
  attr(out, "parameters") <- list(group_by = group_by, split_by = split_by,
                                  within = within, min_dp_site = min_dp_site,
                                  min_cells = min_cells, common_sites = common_sites,
                                  weighted = weighted, n_sites_input = length(site_idx))
  attr(out, "min_attainable_p_paired") <- .min_attainable_p_paired(n_pairs)

  if (verbose) {
    cat(sprintf("Returned %d pseudobulks across %d %s units.\n",
                nrow(out), n_pairs, split_by))
    cat(sprintf("Paired test over %d units floors at p = %.4g.\n",
                n_pairs, .min_attainable_p_paired(n_pairs)))
  }
  out
})


#' @title findDEAlleleFraction: Per-Site Alt-Allele Fraction Between Groups
#' @name findDEAlleleFraction
#'
#' @description
#' Tests, at each site, whether the fraction of reads carrying the alternative
#' base differs between two groups. This is the statistic for RNA-level variation
#' - editing rate, allele-specific expression, heteroplasmy - as distinct from
#' \code{findSNPsByGroup()} (genotype) and \code{findDESNPs()} (read depth).
#'
#' Cells are aggregated to one pseudobulk per `split_by` unit per arm before
#' testing, so the unit of observation is the sample rather than the cell.
#'
#' @section Statistical ceiling:
#' The function reports the smallest p-value attainable given the number of
#' units, and warns when that ceiling cannot clear Bonferroni over the number of
#' sites tested. This is a permutation-count limit that no amount of depth or
#' cell number can overcome; when it binds, the correct move is either to
#' restrict to a small candidate site list or to use
#' \code{computeAlleleFractionIndex()}, which spends only one test.
#'
#' @param ident.1 Character. First group.
#' @param ident.2 Character or NULL. Second group; NULL means all other cells.
#' @param group_by Character. Identity column. Defaults to current identity.
#' @param split_by Character. Unit of observation, default "sample_id".
#' @param sites Optional site restriction (indices, rs IDs, or "chr:pos").
#' @param min_dp Numeric. Minimum pseudobulk depth per arm per unit.
#' @param min_units Numeric. Minimum units per arm for a site to be tested.
#' @param paired Logical. If TRUE, units contributing both arms are tested with
#'   a paired signed-rank test on the per-unit fraction difference.
#' @param test Character. "wilcoxon" (default, distribution-free) or
#'   "betabinom" (pooled likelihood ratio with method-of-moments overdispersion).
#' @param p_adjust Character. Passed to \code{p.adjust}.
#' @param verbose Logical. Print progress.
#'
#' @return A data frame of per-site results ordered by p-value, with attributes
#'   carrying the attainable ceiling and the Bonferroni threshold.
#' @export
variantCell$set("public", "findDEAlleleFraction", function(ident.1,
                                                           ident.2 = NULL,
                                                           group_by = NULL,
                                                           split_by = "sample_id",
                                                           sites = NULL,
                                                           min_dp = 10,
                                                           min_units = 3,
                                                           paired = TRUE,
                                                           test = c("wilcoxon", "betabinom"),
                                                           p_adjust = "BH",
                                                           verbose = TRUE) {
  test <- match.arg(test)
  if (is.null(self$snp_database)) stop("No SNP database. Run buildSNPDatabase() first.", call. = FALSE)
  meta <- self$snp_database$cell_metadata
  if (is.null(group_by)) group_by <- self$current_project_ident
  if (is.null(group_by)) stop("No identity set. Pass group_by or run setProjectIdentity().", call. = FALSE)
  if (!group_by %in% colnames(meta)) stop(sprintf("Column '%s' not found", group_by), call. = FALSE)
  if (!split_by %in% colnames(meta)) stop(sprintf("Column '%s' not found", split_by), call. = FALSE)

  arm <- ifelse(meta[[group_by]] == ident.1, "g1",
                if (is.null(ident.2)) "g2" else ifelse(meta[[group_by]] == ident.2, "g2", NA))
  if (!any(arm == "g1", na.rm = TRUE)) stop(sprintf("No cells in group '%s'", ident.1), call. = FALSE)
  if (!any(arm == "g2", na.rm = TRUE)) stop("Second group is empty", call. = FALSE)

  site_idx <- .resolve_sites(self$snp_database, sites)
  AD <- self$snp_database$ad_matrix[site_idx, , drop = FALSE]
  DP <- self$snp_database$dp_matrix[site_idx, , drop = FALSE]

  key <- paste(meta[[split_by]], arm, sep = "\r")
  key[is.na(arm) | is.na(meta[[split_by]])] <- NA
  if (verbose) cat(sprintf("Aggregating %d sites to pseudobulk...\n", length(site_idx)))
  agg <- .aggregate_pseudobulk(AD, DP, key)

  fields <- do.call(rbind, strsplit(colnames(agg$DP), "\r", fixed = TRUE))
  unit <- fields[, 1]; which_arm <- fields[, 2]
  u1 <- which(which_arm == "g1"); u2 <- which(which_arm == "g2")
  shared <- intersect(unit[u1], unit[u2])

  if (paired && length(shared) < 2) {
    warning("Fewer than 2 units contribute both arms; falling back to unpaired.", call. = FALSE)
    paired <- FALSE
  }

  if (paired) {
    c1 <- u1[match(shared, unit[u1])]
    c2 <- u2[match(shared, unit[u2])]
    n_eff <- length(shared)
    ceiling_p <- .min_attainable_p_paired(n_eff)
  } else {
    c1 <- u1; c2 <- u2
    n_eff <- min(length(c1), length(c2))
    ceiling_p <- .min_attainable_p(length(c1), length(c2))
  }

  ok1 <- agg$DP[, c1, drop = FALSE] >= min_dp
  ok2 <- agg$DP[, c2, drop = FALSE] >= min_dp
  testable <- if (paired) rowSums(ok1 & ok2) >= min_units
              else (rowSums(ok1) >= min_units & rowSums(ok2) >= min_units)

  rows <- which(testable)
  if (length(rows) == 0) {
    warning("No sites met the coverage criteria.", call. = FALSE)
    return(data.frame())
  }
  if (verbose) cat(sprintf("Testing %d sites (%s, %s)...\n", length(rows),
                           test, if (paired) "paired" else "unpaired"))

  bonf <- 0.05 / length(rows)
  if (!is.na(ceiling_p) && ceiling_p > bonf) {
    warning(sprintf(paste0("Best attainable p (%.4g) cannot clear Bonferroni (%.3g) over %d sites - ",
                           "short by %.0fx. Per-site significance is unreachable at this n. ",
                           "Restrict `sites` to a candidate list, or use computeAlleleFractionIndex()."),
                    ceiling_p, bonf, length(rows), ceiling_p / bonf), call. = FALSE)
  }

  rho <- if (test == "betabinom") .estimate_rho(agg$AD[rows, c(c1, c2)], agg$DP[rows, c(c1, c2)]) else NA_real_

  res <- lapply(rows, function(i) {
    a1 <- agg$AD[i, c1]; d1 <- agg$DP[i, c1]
    a2 <- agg$AD[i, c2]; d2 <- agg$DP[i, c2]
    if (paired) {
      keep <- d1 >= min_dp & d2 >= min_dp
      a1 <- a1[keep]; d1 <- d1[keep]; a2 <- a2[keep]; d2 <- d2[keep]
    } else {
      k1 <- d1 >= min_dp; k2 <- d2 >= min_dp
      a1 <- a1[k1]; d1 <- d1[k1]; a2 <- a2[k2]; d2 <- d2[k2]
    }
    f1 <- a1 / d1; f2 <- a2 / d2
    p <- if (test == "betabinom") {
      .betabinom_lrt(a1, d1, a2, d2, rho)
    } else if (paired) {
      if (length(f1) < 2 || all(f1 == f2)) NA_real_
      else suppressWarnings(stats::wilcox.test(f1, f2, paired = TRUE)$p.value)
    } else {
      if (length(f1) < 1 || length(f2) < 1) NA_real_
      else suppressWarnings(stats::wilcox.test(f1, f2)$p.value)
    }
    c(mean(f1), mean(f2), sum(d1), sum(d2), length(f1), p)
  })
  res <- do.call(rbind, res)

  keep_rows <- site_idx[rows]
  info <- self$snp_database$snp_info[keep_rows, , drop = FALSE]
  # gene_name is carried by snp_metrics / snp_annotations rather than snp_info.
  gene <- NA_character_
  for (tb in c("snp_metrics", "snp_annotations")) {
    tbl <- self$snp_database[[tb]]
    if (!is.null(tbl) && "gene_name" %in% colnames(tbl) &&
        nrow(tbl) == nrow(self$snp_database$snp_info)) {
      gene <- tbl$gene_name[keep_rows]
      break
    }
  }
  out <- data.frame(
    CHROM      = info$CHROM,
    POS        = info$POS,
    rs_id      = if (!is.null(info$rs_id)) info$rs_id else NA_character_,
    gene_name  = gene,
    af_1       = res[, 1],
    af_2       = res[, 2],
    af_diff    = res[, 1] - res[, 2],
    dp_1       = res[, 3],
    dp_2       = res[, 4],
    n_units    = res[, 5],
    p_value    = res[, 6],
    stringsAsFactors = FALSE
  )
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust)
  out <- out[order(out$p_value, -abs(out$af_diff)), , drop = FALSE]
  rownames(out) <- NULL

  attr(out, "min_attainable_p")    <- ceiling_p
  attr(out, "bonferroni_threshold") <- bonf
  attr(out, "significance_reachable") <- !is.na(ceiling_p) && ceiling_p <= bonf
  attr(out, "parameters") <- list(ident.1 = ident.1, ident.2 = ident.2,
                                  group_by = group_by, split_by = split_by,
                                  paired = paired, test = test, rho = rho,
                                  min_dp = min_dp, min_units = min_units,
                                  n_units = n_eff)
  if (verbose) {
    cat(sprintf("Best attainable p at n=%d: %.4g   Bonferroni: %.3g   reachable: %s\n",
                n_eff, ceiling_p, bonf,
                if (!is.na(ceiling_p) && ceiling_p <= bonf) "YES" else "NO"))
  }
  out
})
