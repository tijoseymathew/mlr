#' @title Visualize a prediction object. Right now supports only PHMRegr tasks
#'
#' @description
#' Plots the true value and the response in prediction via \code{\link[ggplot2]{ggplot}}.
#'
#' FIXME: Implement for general classification & regression problems. Extend from plotLearnerPrediction.R
#'
#' For phm regression, the true value and the predicted response of target is plotted for each seq.id in separate facets.
#' 
#' The plot title displays the learner name, its parameters, the training performance and the cross-validation performance.
#'
#' @param prediction [\code{Prediction}]\cr
#'   The prediction object to be plotted
#' @param seq.ids [\code{character(1)} or \code{numeric(1)}]
#'   For phm regression: Either a vector of seq.id to be plotted or a numeric fraction of seq.id to be randomly sampled
#'   Default is 1.
#' @param linesize [\code{numeric(1)}]\cr
#'   Linesize for ggplot2 \code{\link[ggplot2]{geom_point}} for data points.
#'   Default is 1.
#' @return The ggplot2 object.
#' @export
plotPrediction = function(prediction, linesize = 1, greyscale = FALSE, seq.ids = 1) {

  assertClass(prediction, "Prediction")
  if (!checkClass(prediction, "PredictionPHMRegr"))
    stop("Only PHM Regr plotting supported for now")
  td = prediction$task.desc
  sid = td$seq.id; tid = td$order.by; tar = td$target

  pltDF = prediction$data
  pltDF = melt(pltDF, measure.vars = c("truth", "response"), 
               variable.name = "Type")
  if (is.numeric(seq.ids)) {
    assertNumeric(seq.ids, lower = 0, upper = 1)
    u_sids = unique(pltDF[[sid]])
    seq.ids = sample(u_sids, round(length(u_sids) * seq.ids))
  }
  pltDF = pltDF[pltDF[[sid]] %in% seq.ids, ]
  p = ggplot(pltDF, aes_string(x = tid, y = "value"))+
    xlab(tid)+ylab(tar)+
    geom_line(aes(col = Type), size = linesize)+
    facet_grid(as.formula(paste(sid, "~.")))

  # set title
  title = sprintf("Predictions on %s", td$id)
  if (is.character(seq.ids))
    title = paste("Partial", title)
  p = p + ggtitle(title)

  return(p)
}