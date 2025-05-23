
#' @export
variantCell <- R6Class("variantCell",
                       public = list(
                         samples = list(),
                         snp_database = data.frame(),
                         current_project_ident = NULL,

                         initialize = function() {
                           self$samples <- list()
                           self$snp_database <- list()
                           self$current_project_ident <- NULL
                         }))
