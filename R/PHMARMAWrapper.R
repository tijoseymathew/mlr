#' @title PHM with ARMA features
#'
#' @description
#' Wraps a PHM learner that trains and predicts on ARMA features
#'
#' @template arg_learner
#' @param lag.ar [\code{integer(1)}]\cr
#'   Auto-regressive (target) lag to be incorporated
#'   Default in 0
#' @param lag.ma [\code{integer(1)}]\cr
#'   Moving Average (feature) lag to be incorporated
#'   Default is 0
#' @param fill.arma [\code{character(1)}]\cr
#'   Value to be filled at begining of lagged columns
#'   \code{"zero"} for 0 and \code{"mean"} for mean value
#' @param ar.stage [\code{character(1)}]\cr
#'   The stage of using auto-regressive features
#'   \code{"single-step"} for using auto-regressive \code{target} values
#'   \code{"two-step"} for using auto-regressive features of \code{target} predicted in first step using only \code{features}
#' @template ret_learner
#' @export
#' @family wrapper
#' @examples
#' lrn = makeLearner("regr.ksvm")
#' lrn = makePHMARMAWrapper(lrn)
#' mod = train(lrn, phm.task)
#' predictions = predict(mod, newdata = getTaskData(phm.task))
makePHMARMAWrapper = function(learner) {
  checkLearner(learner, "phmregr")
  ps = makeParamSet(
    makeIntegerLearnerParam(id = "lag.ar", lower = 0L, default = 0L),
    makeIntegerLearnerParam(id = "lag.ma", lower = 0L, default = 0L),
    makeDiscreteLearnerParam(id = "fill.arma", values = c("mean", "zero"), default = "mean"),
    makeDiscreteLearnerParam(id = "ar.stage", values = c("single-step", "two-step"), default = "single-step")
  )
  lrn = makeBaseWrapper(id = paste(learner$id, "phmarma", sep = "."),
                        type = "phmregr", next.learner = learner,
                        par.set = ps,
                        learner.subclass = "PHMARMAWrapper",
                        model.subclass = "PHMARMAModel")
  lrn$properties = learner$properties
  return(lrn)
}

#' @export
trainLearner.PHMARMAWrapper = function(.learner, .task, .subset, .weights = NULL, 
                                       lag.ar = 0L, lag.ma = 0L, fill.arma = "mean",
                                       ar.stage = "single-step", 
                                       ...) {
  assertClass(.task, "PHMRegrTask")
  td = .task$task.desc
  args = list(seq.id = td$seq.id, order.by = td$order.by, target = td$target,
              lag.ma = lag.ma, fill.arma = fill.arma)
  if (ar.stage == "single-step") {
    args$lag.ar = lag.ar
    arma = makeARMADF(getTaskData(.task, .subset), args)
    
    task = makePHMRegrTask(id = paste0(getTaskId(.task), ".PHMARMA"),
                           data = arma$data, target = td$target,
                           seq.id = td$seq.id, order.by = td$order.by)
    
    n.mdl = train(.learner$next.learner, task, weights = .weights) 
    arma.args = arma$args
  } else if (ar.stage == "two-step") {
    tskData = getTaskData(.task, .subset)
    n.mdl = list()
    arma.args = list()
    
    args$lag.ar = 0L
    arma = makeARMADF(tskData, args)
    task = makePHMRegrTask(id = paste0(getTaskId(.task), ".PHMARMA_Step1"),
                           data = arma$data, target = td$target,
                           seq.id = td$seq.id, order.by = td$order.by)
    n.mdl[[1]] = train(.learner$next.learner, task, weights = .weights)
    arma.args[[1]] = arma$args
    # Stage-2
    y = tskData[[td$target]]
    tskData[[td$target]] = predict(n.mdl[[1]], task)$data$response
    # FIXME: By def arma ar will not use value at t. Not logical here
    args$lag.ar = lag.ar
    arma = makeARMADF(tskData, args)
    arma$data[[td$target]] = y 
    task = makePHMRegrTask(id = paste0(getTaskId(.task), ".PHMARMA_Step1"),
                           data = arma$data, target = td$target,
                           seq.id = td$seq.id, order.by = td$order.by)
    n.mdl[[2]] = train(.learner$next.learner, task, weights = .weights)
    arma.args[[2]] = arma$args
  }
  
  m = makeWrappedModel.Learner(learner = .learner, 
                               learner.model = list(regr.model = n.mdl, arma.args = arma.args, ar.stage = ar.stage),
                               task.desc = getTaskDesc(.task),
                               subset = .subset,
                               features = getTaskFeatureNames(.task),
                               factor.levels = mlr:::getTaskFactorLevels(.task),
                               time = 0)
  makeChainModel(next.model = m, 
                 cl = "PHMViaRegressionModel")
}

#' @export
predictLearner.PHMARMAWrapper = function(.learner, .model, .newdata, ...) {
  args = removeFromDots(names(.learner$par.vals), ...)
  
  mdl = .model$learner.model$next.model$learner.model
  arma.args = mdl$arma.args
  if (mdl$ar.stage == "two-step") {
    data = makeARMADF(.newdata, arma.args[[1]])$data
    ret = do.call( predictLearner, 
                   c( list(.learner = .learner$next.learner, 
                           .model = mdl$regr.model[[1]], 
                           .newdata = data), 
                      args) )
    .newdata[[arma.args[[1]]$target]] = ret$y
    data = makeARMADF(.newdata, arma.args[[2]])$data
    data[[arma.args[[1]]$target]] = NULL
    ret = do.call( predictLearner, 
                   c( list(.learner = .learner$next.learner, 
                           .model = mdl$regr.model[[2]], 
                           .newdata = data), 
                      args) )
    
  } else if (arma.args$lag.ar == 0 ) {
    data = makeARMADF(.newdata, arma.args)$data
    ret = do.call(predictLearner, 
                  c( list(.learner = .learner$next.learner, 
                          .model = mdl$regr.model, 
                          .newdata = data), 
                     args)    
    )  
  } 
  else {
    sid = arma.args$seq.id; tid = arma.args$order.by; tar = arma.args$target
    ret = as.data.table(.newdata[, c(sid, tid)])
    ret$y = numeric()
    
    # Get MA data
    ma.args = mdl$arma.args
    ma.args$lag.ar = 0
    dat = as.data.table(makeARMADF(.newdata, ma.args)$data)
    
    dat[order(get(tid)), .time := 1:.N, by = list(get(sid))]
    # AR initialization
    cns = paste0(tar, "_lag_", 1:arma.args$lag.ar)
    seqs = unique(dat[[sid]])
    arVals = data.frame( matrix(1, nrow=length(seqs)) %*% t(unlist(arma.args[cns])) )
    arVals[[sid]] = seqs
    # One-step ahead prediction
    for (step in 1:max(dat$.time)) {
      pdat = as.data.frame(dat[.time == step])
      pdat = merge(pdat, arVals, by = sid, all = FALSE)
      pdat$.time = NULL
      y = as.data.table( do.call(predictLearner, 
                                 c( list(.learner = .learner$next.learner, 
                                         .model = mdl$regr.model, 
                                         .newdata = pdat), 
                                    args)) )
      ret[y, ':=' (y = i.y), on = c(sid, tid)]
      # Shift initial values
      arVals = arVals[arVals[[sid]] %in% y[[sid]], ]
      arVals = cbind(y$y, arVals[, cns, drop = FALSE][, -length(cns), drop=FALSE])
      colnames(arVals) = cns
      arVals[[sid]] = y[[sid]]
    }
  }
  as.data.frame(ret)
}

# Takes in a data frame of seq.id, order.by, target, and args specifying the parameters and returns a ARMA version
makeARMADF = function(data, args) {
  sid = args$seq.id; tid = args$order.by; tar = args$target
  
  maCols = setdiff(colnames(data), c(sid, tid, tar))
  l.maCns = unlist(lapply(maCols, function(cn) paste0(cn, "_lag_", 0:args$lag.ma)))
  
  data = as.data.table(data)
  if (args$lag.ar > 0) { # Present value should not be used as feature
    l.arCns = paste0(tar, "_lag_", 1:args$lag.ar)
    data[order(get(tid)),
         (l.arCns) := shift(.SD, n=1:args$lag.ar),
         by = get(sid), .SDcols = tar]
  }
  data[order(get(tid)),
       (l.maCns) := shift(.SD, n=0:args$lag.ma),
       by = get(sid), .SDcols = maCols]
  
  fillFn = switch(args$fill.arma,
                  "mean" = function(x) mean(x, na.rm = TRUE),
                  "zero" = function(x) 0)
  lcns = grep("_lag_", names(data), value = TRUE)
  data = data[, intersect(c(lcns, sid, tid, tar), colnames(data)), 
              with = FALSE]
  for (col in lcns) {
    if (col %in% names(args)) {
      fill = args[[col]]
    } else {
      fill = fillFn(data[[col]])
      args[[col]] = fill
    }
    data[[col]][is.na(data[[col]])] = fill
  }
  list(data = as.data.frame(data), args = args)
}
