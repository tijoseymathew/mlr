#' @title Plot a benchmark summary.
#'
#' @description
#' Creates a facetted plot for each where the true response and predicted response from the learners are plotted out
#'
#' @template arg_bmr
#' @param prediction [\code{Prediction}]\cr
#'   The prediction object to be plotted
#' @param seq.ids [\code{character(1)} or \code{numeric(1)}]
#'   For phm regression: Either a vector of seq.id to be plotted or a numeric fraction of seq.id to be randomly sampled
#'   Default is 1.
#' @param linesize [\code{numeric(1)}]\cr
#'   Linesize for ggplot2 \code{\link[ggplot2]{geom_point}} for data points.
#'   Default is 1.
#' @family benchmark
#' @family plot
#' @export
#' @examples
#' # see benchmark
plotBMRPrediction = function(bmr, linesize = 1, seq.ids = 1) {
  assertClass(bmr, "BenchmarkResult")
  
  if (length(getBMRTaskIds(bmr)) > 1)
    stop("Only one task supported for now!")
  # Test to see if phmregr
  tds = getBMRTaskDescs(bmr)
  tds = unlist(tds)
  if (!all(tds[grepl("type", names(tds))] == "phmregr"))
    stop("Only PHM Regr plotting supported for now!")
  
  sid = unique(tds[grepl("seq.id", names(tds))])
  tid = unique(tds[grepl("order.by", names(tds))])
  tar = unique(tds[grepl("target", names(tds))])
  
  pltDF = setDT(getBMRPredictions(bmr, as.df = TRUE))
  trDT = unique(pltDF[, c("task.id", "id", sid, tid, "truth", "iter", "set"),
                      with = FALSE])
  trDT$learner.id = "truth"
  pltDF$truth = NULL
  names(pltDF)[names(pltDF) == "response"] = "value"
  names(trDT)[names(trDT) == "truth"] = "value"
  pltDF = rbind(pltDF, trDT)
  
  if (is.numeric(seq.ids)) {
    assertNumeric(seq.ids, lower = 0, upper = 1)
    u_sids = unique(pltDF[[sid]])
    seq.ids = sample(u_sids, round(length(u_sids) * seq.ids))
  }
  pltDF = pltDF[pltDF[[sid]] %in% seq.ids]
  p = ggplot(pltDF, aes_string(x = tid, y = "value"))+
    xlab(tid)+ylab(tar)+
    geom_line(aes(col = learner.id), size = linesize)+
    facet_grid(as.formula(paste(sid, "~.")))
  
  return(p)
}