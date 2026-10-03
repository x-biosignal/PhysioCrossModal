#' The multimodal container
#'
#' \pkg{PhysioCrossModal} used to define its own \code{MultiPhysioExperiment}
#' class, holding \code{experiments} / \code{alignment} / \code{sampleMap} /
#' \code{couplingResults}. That class no longer exists here: the ecosystem has a
#' single multimodal container, defined in \pkg{PhysioExperiment}, and this
#' package is an analysis package that takes it as input.
#'
#' The constructor and every accessor keep their names and return types, so
#' existing code runs unchanged:
#' \code{MultiPhysioExperiment()}, \code{experiments()}, \code{modalities()},
#' \code{samplingRates()}, \code{nModalities()}, \code{alignment()},
#' \code{alignment<-}, \code{length()}, \code{names()}, \code{[}, \code{[[}.
#' \code{alignment()} is now a view computed from the container's shared clock
#' rather than a separately stored table, so the alignment and the signals
#' cannot drift apart.
#'
#' Containers saved under the old class are converted by
#' \code{\link[PhysioExperiment]{migrateContainer}} or read with
#' \code{\link[PhysioExperiment]{readPhysioRDS}}.
#'
#' @name MultiPhysioExperiment-crossmodal
#' @seealso \code{\link[PhysioExperiment]{MultiPhysioExperiment}},
#'   \code{\link{couplingAnalysis}}, \code{\link{couplingResults}}
NULL
