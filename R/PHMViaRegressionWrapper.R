#' @title PHM via regression wrapper
#'
#' @description
#' Wraps a regression learner to handle PHM tasks by removing the seq.id and order.by columns from the data. 
#' Prediction will re-assign original seq.id and order.by
#'
#' @template arg_learner
#' @template ret_learner
#' @export
#' @family wrapper
#' @examples
#' lrn = makeLearner("regr.ksvm")
#' lrn = makePHMViaRegressionWrapper(lrn)
#' mod = train(lrn, phm.task)
#' predictions = predict(mod, newdata = getTaskData(phm.task))
makePHMViaRegressionWrapper = function(learner) {
  checkLearner(learner, "regr")
  lrn = makeBaseWrapper(id = paste(learner$id, "phmviaregr", sep = "."),
                        type = "phmregr", next.learner = learner, 
                        learner.subclass = "PHMViaRegressionWrapper",
                        model.subclass = "PHMViaRegressionModel")
  lrn$properties = learner$properties
  return(lrn)
}

#' @export
trainLearner.PHMViaRegressionWrapper = function(.learner, .task, .subset, .weights = NULL, ...) {
  assertClass(.task, "PHMRegrTask")
  td = .task$task.desc
  task = dropFeatures(subsetTask(.task, .subset), c(td$seq.id, td$order.by))
  task = makeRegrTask(id = paste0(getTaskId(.task), ".PHMViaRegression"),
                      data = getTaskData(task),
                      target = getTaskTargetNames(task))
  n.mdl = train(.learner$next.learner, task, weights = .weights)
  m = makeWrappedModel.Learner(learner = .learner, 
                               learner.model = list(regr.model = n.mdl, 
                                                    seq.id = td$seq.id, 
                                                    order.by = td$order.by),
                               task.desc = getTaskDesc(.task),
                               subset = .subset,
                               features = getTaskFeatureNames(.task),
                               factor.levels = mlr:::getTaskFactorLevels(.task),
                               time = 0)
  makeChainModel(next.model = m, 
                 cl = "PHMViaRegressionModel")
}

#' @export
predictLearner.PHMViaRegressionWrapper = function(.learner, .model, .newdata, ...) {
  args = removeFromDots(names(.learner$par.vals), ...)
  
  mdl = .model$learner.model$next.model$learner.model
  y = do.call(predictLearner, 
              c( list(.learner = .learner$next.learner, 
                      .model = mdl$regr.model, 
                      .newdata = .newdata[, setdiff(colnames(.newdata), 
                                                    c(mdl$seq.id, mdl$order.by)),
                                          drop = FALSE]), # Converts to numeric if only one feature is present. Hence drop = FALSE
                 args)
  )
  ret = list()
  ret[[mdl$seq.id]] = .newdata[[mdl$seq.id]]
  ret[[mdl$order.by]] = .newdata[[mdl$order.by]]
  ret[["y"]] = y
  data.frame(ret)
}
