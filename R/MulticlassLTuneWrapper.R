# MulticlassWrapper with LearnerTuning
#' @export
makeMulticlassLTuneWrapper = function(learner, mcw.method = "onevsrest", tune.params = list()) {
  learner = checkLearner(learner)
  ps = makeParamSet(
    makeUntypedLearnerParam(id = "mcw.method", default = "onevsrest"),
    makeUntypedLearnerParam(id = "tune.params", default = list())
  )
  assert(
    checkChoice(mcw.method, c("onevsrest", "onevsone")),
    checkFunction(mcw.method, args = "task")
  )
  pv = list(mcw.method = mcw.method, tune.params = tune.params)
  id = stri_paste(learner$id, "multiclassLTune", sep = ".")
  
  x = makeHomogeneousEnsemble(id = id, type = "classif", next.learner = learner,
                              package = learner$package,  par.set = ps, par.vals = pv,
                              learner.subclass = "MulticlassLTuneWrapper", model.subclass = "MulticlassLTuneModel")
  x = setPredictType(x, predict.type = "response")
  return(x)
}

#' @export
trainLearner.MulticlassLTuneWrapper = function(.learner, .task, .subset = NULL, .weights = NULL, mcw.method, tune.params, ...) {
  .task = subsetTask(.task, .subset)
  y = getTaskTargets(.task)
  cm = buildCMatrix(mcw.method, .task)
  x = multi.to.binary(y, cm)
  tLearner = do.call(makeTuneWrapper, c(list(learner=.learner$next.learner), tune.params))
  tLearner = list(next.learner = tLearner) # Hack to use doMulticlassTrainIteration
  args = list(x = x, learner = tLearner, task = .task, weights = .weights)
  parallelLibrary("mlr", master = FALSE, level = "mlr.ensemble", show.info = FALSE)
  exportMlrOptions(level = "mlr.ensemble")
  models = parallelMap(i = seq_along(x$row.inds), doMulticlassTrainIteration,
                       more.args = args, level = "mlr.ensemble")
  m = makeHomChainModel(tLearner, models)
  m$cm = cm
  return(m)
}

predictLearner.MulticlassLTuneWrapper = function(.learner, .model, .newdata, ...) {
  predictLearner.MulticlassWrapper(.learner, .model, .newdata, ...)
}

#' @export
getLearnerProperties.MulticlassLTuneWrapper = function(learner){
  getLearnerProperties.MulticlassWrapper(learner)
}