#' @title Extend learner as an ensemble trained on subset of seq.ids.
#'
#' @description
#' Learner is trained on subsets of sequences for the PHM task to form an ensemble of models
#' Models can easily be accessed via \code{\link{getLearnerModel}}.
#'
#' Models are trained either on each seq.id individually or a group of seqs
#'
#' For prediction the mean value and the standard deviations across predictions is computed.
#'
#' Note that the passed base learner must always have \code{predict.type = 'response'},
#' while the PHMEnsembleWrapper can estimate probabilities and standard errors, so it can
#' be set, e.g., to \code{predict.type = 'prob'}. For this reason, when you call
#' \code{\link{setPredictType}}, the type is only set for the PHMEnsembleWrapper, not passed
#' down to the inner learner.
#'
#' FIXME: Add option to group sequences based on length?
#' @template arg_learner
#' @template ret_learner
#' @family wrapper
#' @export
makePHMEnsembleWrapper = function(learner) {
  learner = checkLearner(learner, type = c("phmregr"))
  id = stri_paste(learner$id, "ensembled", sep = ".")
  packs = learner$package
  makeHomogeneousEnsemble(id, learner$type, learner, packs, learner.subclass = "PHMEnsembleWrapper", model.subclass = "PHMEnsembleModel")
}

#' @export
print.PHMEnsembleModel = function(x, ...) {
  s = capture.output(print.WrappedModel(x))
  u = sprintf("PHM Ensemble Learner: %s", class(x$learner$next.learner)[1L])
  s = append(s, u, 1L)
  lapply(s, catf)
}

#' @export
trainLearner.PHMEnsembleWrapper = function(.learner, .task, .subset, .weights = NULL, ...) {
  
  assertClass(.task, "PHMRegrTask")
  .task = subsetTask(.task, subset = .subset)
  seq.ids = unique( 
    getTaskData(.task, features = getTaskDesc(.task)$seq.id, target.extra = TRUE)$data[, 1]
  )
  
  args = list(task = .task, learner = .learner, weights = .weights)
  parallelLibrary("mlr", master = FALSE, level = "mlr.ensemble", show.info = FALSE)
  exportMlrOptions(level = "mlr.ensemble")
  models = parallelMap(doPHMEnsembleTrainIteration, seq.ids, more.args = args, level = "mlr.ensemble")
  makeHomChainModel(.learner, models)
}

doPHMEnsembleTrainIteration = function(seq.id, task, learner, weights) {
  setSlaveOptions()
  sids = getTaskData(task, features = getTaskDesc(task)$seq.id, target.extra = TRUE)$data[, 1]
  ix = (sids == seq.id)
  train(learner$next.learner, task, subset = ix, weights = weights[ix])
}

#' @export
predictLearner.PHMEnsembleWrapper = function(.learner, .model, .newdata, ...) {
  models = getLearnerModel(.model, more.unwrap = FALSE)
  p = asMatrixCols(lapply(models, function(m) {
    predict(m, newdata = .newdata, ...)$data$response
  }))
  ret = predict(models[[1]], newdata = .newdata, ...)$data
  ret$response = NULL
  ret$y = rowMeans(p)
  ret
}
