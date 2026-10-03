# The multimodal container is defined once, in PhysioExperiment, which this
# package depends on directly. These re-exports keep `PhysioCrossModal::`
# qualified calls working for code written when the container lived here.
# Nothing is redefined: each name below is the same object.

#' @importFrom PhysioExperiment MultiPhysioExperiment
#' @export
PhysioExperiment::MultiPhysioExperiment

#' @importFrom PhysioExperiment experiments
#' @export
PhysioExperiment::experiments

#' @importFrom PhysioExperiment `experiments<-`
#' @export
PhysioExperiment::`experiments<-`

#' @importFrom PhysioExperiment modalities
#' @export
PhysioExperiment::modalities

#' @importFrom PhysioExperiment samplingRates
#' @export
PhysioExperiment::samplingRates

#' @importFrom PhysioExperiment nModalities
#' @export
PhysioExperiment::nModalities

#' @importFrom PhysioExperiment alignment
#' @export
PhysioExperiment::alignment

#' @importFrom PhysioExperiment `alignment<-`
#' @export
PhysioExperiment::`alignment<-`

#' @importFrom PhysioExperiment streams
#' @export
PhysioExperiment::streams

#' @importFrom PhysioExperiment commonClock
#' @export
PhysioExperiment::commonClock

#' @importFrom PhysioExperiment timeWindow
#' @export
PhysioExperiment::timeWindow

#' @importClassesFrom PhysioExperiment MultiPhysioExperiment
NULL
