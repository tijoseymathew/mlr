#' @title PHM with Health Index Transformation
#'
#' @description
#' Wraps a PHM learner that trains to target a given health index that is scaled back during prediction
#'
#' @template arg_learner
#' @param hiFns [\code{list(2)}]\cr
#'   List with trainFn and predictFn for HI model.
#'   trainFn: Function that returns a model and health index for each observation in task
#'   predictFn: Function that computes the RUL for data predicted using arg_learner
#' @template ret_learner
#' @export
#' @family wrapper
#' @examples
#' lrn = makePHMViaRegressionWrapper(makeLearner("regr.ksvm"))
#' lrn = makePHMHIWrapper(lrn)
#' mod = train(lrn, phm.task)
#' predictions = predict(mod, newdata = getTaskData(phm.task))
makePHMHIWrapper = function(learner) {
  checkLearner(learner, "phmregr")
  ps = makeParamSet(
    makeUntypedLearnerParam("hiFn", tunable = FALSE)
  )
  lrn = makeBaseWrapper(id = paste(learner$id, "hiFn", sep = "."),
                        type = "phmregr", next.learner = learner,
                        par.set = ps, 
                        learner.subclass = "PHMHIWrapper",
                        model.subclass = "PHMHIModel")
  lrn$properties = learner$properties
  return(lrn)
}

#' @export
trainLearner.PHMHIWrapper = function(.learner, .task, .subset, .weights = NULL, 
                                     hiFns = linearHI(), ...) {
  assertClass(.task, "PHMRegrTask")
  assertFunction(hiFns$trainFn, args = "task")
  assertFunction(hiFns$predictFn, args = c("model", "data"))
  
  task = subsetTask(.task, .subset)
  hi.op = hiFns$trainFn(task)
  assertNumeric(hi.op$hi, len = getTaskSize(task), any.missing = FALSE)
  
  data = getTaskData(task, target.extra = TRUE)$data
  data$.hi = hi.op$hi
  td = getTaskDesc(task)
  task = makePHMRegrTask(id = paste0(getTaskId(.task), ".hi"),
                         data = data, target = ".hi",
                         seq.id = td$seq.id, order.by = td$order.by)
  next.mdl = train(.learner$next.learner, task, weights = .weights)
  
  m = makeWrappedModel.Learner(learner = .learner, 
                               learner.model = list(next.model = next.mdl, 
                                                    hi.model = hi.op$model, predictFn = hiFns$predictFn),
                               task.desc = getTaskDesc(.task),
                               subset = .subset,
                               features = getTaskFeatureNames(.task),
                               factor.levels = mlr:::getTaskFactorLevels(.task),
                               time = 0)
  ret = makeChainModel(next.model = m, cl = "PHMHIModel")
}

#' @export
predictLearner.PHMHIWrapper = function(.learner, .model, .newdata, ...) {
  args = removeFromDots(names(.learner$par.vals), ...)
  
  mdl = .model$learner.model$next.model$learner.model
  ret = do.call( predictLearner, 
                 c( list(.learner = .learner$next.learner, 
                         .model = mdl$next.model,
                         .newdata = .newdata), 
                    args) )
  y = mdl$predictFn(mdl$hi.model, ret)
  assertNumeric(y, len = nrow(.newdata), any.missing = FALSE)
  ret$y = y 
  ret
}


# HI Functions ------------------------------------------------------------

#' @title HI model functions for PHMHIWrapper
#'
#' @description
#' HI from features on a window at start and end of each seq.id
#' 
#' @param type [\code{character(1)}]\cr
#'   Either "linear" or "exponential"
#' @param start.n [\code{numeric(1)}]\cr
#'   Number of points at begining to be modelled as HI 1. Default 1
#' @param end.n [\code{numeric(1)}]\cr
#'   Number of points at end to be modelled as HI 0. Default 1
#' @param scaled [\code{logical(1)}]\cr
#'   Should the features be scaled. Default True
#' @param multiplier [\code{numeric(1)}]
#'   Multiplier to be used during prediction. If NULL mean of train RUL will be used. Defaults to NULL
#' @export
linearHI = function(type = "linear", start.n = 1, end.n = 1, scaled = TRUE, multiplier = NULL) {
  linFn = function(trDT, teDT) {
    mdl = lm(hi~., trDT)
    predict(mdl, teDT)
  }
  trainFn = function(task) {
    td = getTaskDesc(task)
    tid = td$order.by; sid = td$seq.id; tar = td$target
    
    dat = setDT( getTaskData(task, target.extra = TRUE)$data )
    if (scaled) {
      cols = setdiff(names(dat), c(tid, sid))
      dat[, (cols) := lapply(.SD, scale), .SDcols = cols]
    }
    
    dat[order(get(tid)), hi := {
      # FIXME: Clip start.n and end.n
      marginDT = rbind(head(.SD, start.n), tail(.SD, end.n))
      marginDT$hi = rep(c(1, 0), times = c(start.n, end.n))
      marginDT[[tid]] = NULL
      
      linFn(marginDT, .SD)
    }, by = sid]
    if (is.null(multiplier))
      multiplier = mean(getTaskTargets(task))
    list(model = multiplier, hi = dat$hi)
  }
  
  predictFn = function(model, data) {
    data$y * model
  }
  list(trainFn = trainFn, predictFn = predictFn)
}