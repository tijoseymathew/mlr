#' @title PHM with stage classification
#'
#' @description
#' Wraps a PHM learner that initially trains a simple classifier to detect normal/abnormal stage and trains PHM learner only on abnormal data
#' @export
makePHMStageClassifierWrapper = function(learner) {
  checkLearner(learner, "phmregr")
  ps = makeParamSet(
    makeUntypedLearnerParam("stage_learner", tunable = FALSE, default = makeLearner("classif.ksvm")),
    makeNumericLearnerParam("target_threshold"),
    makeIntegerLearnerParam("stage_window")
  )
  lrn = makeBaseWrapper(id = paste(learner$id, "stgClass", sep = "."),
                        type = "phmregr", next.learner = learner,
                        par.set = ps, 
                        learner.subclass = "PHMStageClassifierWrapper",
                        model.subclass = "PHMStageClassifierModel")
  lrn$properties = learner$properties
  lrn
}

#' @export
trainLearner.PHMStageClassifierWrapper = function(.learner, .task, .subset, .weights = NULL, 
                                                  stage_learner = makeLearner("classif.ksvm"), 
                                                  target_threshold, # Target below this threshold is considered abnormal
                                                  stage_window, ...) {
  assertClass(.task, "PHMRegrTask")
  checkLearner(stage_learner, "classif")
  assertInteger(stage_window)
  assertNumeric(target_threshold)
  
  task = subsetTask(.task, .subset)
  tid = task$task.desc$order.by; sid = task$task.desc$seq.id; tar = task$task.desc$target
  # Convert to classification task
  dat = getTaskData(task)
  dat$.stage = ifelse(getTaskTargets(task) < target_threshold, ".abnormal", ".normal")
  dat = dat[, !names(dat) %in% c(tid, sid, tar)]
  cl_task = makeClassifTask(data = dat, target = ".stage")
  cl_mdl = train(stage_learner, cl_task)
  # PHM task
  ph_task = subsetTask(task, dat$.stage == ".abnormal")
  next.mdl = train(.learner$next.learner, ph_task, weights = .weights[dat$.stage == ".abnormal"])
  
  m = makeWrappedModel.Learner(learner = .learner, 
                               learner.model = list(next.model = next.mdl, 
                                                    cl_lrn = stage_learner,
                                                    cl_mdl = cl_mdl,
                                                    seq.id = sid,
                                                    order.by = tid,
                                                    stage_window = stage_window,
                                                    target_threshold = target_threshold),
                               task.desc = getTaskDesc(.task),
                               subset = .subset,
                               features = getTaskFeatureNames(.task),
                               factor.levels = mlr:::getTaskFactorLevels(.task),
                               time = 0)
  makeChainModel(next.model = m, cl = "PHMStageClassifierModel")
}

#' @export
predictLearner.PHMStageClassifierWrapper = function(.learner, .model, .newdata, ...) {
  args = removeFromDots(names(.learner$par.vals), ...)
  mdl = .model$learner.model$next.model$learner.model
  
  # Stage classifier
  cl_dat = .newdata[, !names(.newdata) %in% c(mdl$order.by, mdl$seq.id)]
  cl_y = do.call( predictLearner,
                  c( list(.learner = mdl$cl_lrn,
                          .model = mdl$cl_mdl,
                          .newdata = cl_dat),
                     args) )
  # Stage filtering
  dat = data.table(s = .newdata[[mdl$seq.id]],
                   t = .newdata[[mdl$order.by]],
                   stg = (cl_y == ".abnormal"))
  dat[, new_y := {
    d = copy(.SD)
    d$ix = 1:nrow(.SD)
    d = d[order(t)]
    ifelse(cumsum(d$stg) >= mdl$stage_window, ".abnormal", ".normal")[d$ix]
  }, by = s]
  cl_y = dat$new_y
  
  ret = list()
  ret[[mdl$seq.id]] = .newdata[[mdl$seq.id]]
  ret[[mdl$order.by]] = .newdata[[mdl$order.by]]
  ret[["y"]] = rep(mdl$target_threshold, nrow(.newdata))
  ret[["y"]][cl_y == ".abnormal"] = do.call( predictLearner,
                                             c( list(.learner = .learner$next.learner,
                                                     .model = mdl$next.model,
                                                     .newdata = .newdata[cl_y == ".abnormal", , drop=FALSE]),
                                                args) )$y
  data.frame(ret)
}