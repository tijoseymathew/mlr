#' @title PHM with stage classification
#'
#' @description
#' Wraps a PHM learner that initially trains a simple classifier to detect normal/abnormal stage and trains PHM learner only on abnormal data
#' @export
makePHMStageClassifierWrapper = function(learner) {
  checkLearner(learner, "phmregr")
  ps = makeParamSet(
    makeUntypedLearnerParam("stage_learner", tunable = FALSE, default = stagePCAChangePoint)
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
                                                  stage_learner = stagePCAChangePoint, # Returns a model and named list of subsets for each stage. A predict method for model will return subset for each model
                                                  ...) {
  assertClass(.task, "PHMRegrTask")
  
  task = subsetTask(.task, .subset)
  # Stage detection
  stage_obj = stage_learner(task)
  # FIXME: Ensure subsets are consecutive sections of data?
  subs = stage_obj$subsets[!sapply(stage_obj$subsets, is.null)]
  next.mdls = lapply(subs, function(s) train(learner = .learner$next.learner,
                                             task = task,
                                             subset = s,
                                             weights = .weights[s]))
  
  m = makeWrappedModel.Learner(learner = .learner, 
                               learner.model = list(next.model = next.mdls,
                                                    stage.model = stage_obj$model),
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
  
  # Stage detector
  stage_subsets = predict(mdl$stage.model, .newdata)
  if (length(do.call(intersect, unname(stage_subsets))) != 0 | length(do.call(union, unname(stage_subsets))) != nrow(.newdata))
    stopf("Stage function should have all subsets without common elements!")
  ret = lapply(names(mdl$next.model), 
               function(stage) {
                 if (!is.null(stage_subsets[stage]))
                   do.call( predictLearner,
                            c( list(.learner = .learner$next.learner,
                                    .model = mdl$next.model[[stage]],
                                    .newdata = .newdata[stage_subsets[[stage]], ]),
                               args) )
               })
  ret = do.call(rbind, ret)
  ret[order(unlist(stage_subsets)), ]
}

#' @export
stagePCAChangePoint = function(task) {
  td = getTaskDesc(task)
  tid = td$order.by; sid = td$seq.id
  
  dat = as.data.table(getTaskData(task, target.extra = T)$data)
  datMat = as.matrix(dat[, setdiff(names(dat), c(tid, sid)), with=F])
  
  pca_mdl = prcomp(datMat, center = T, scale. = T)
  
  mdl = makeS3Obj(classes = "stagePCACP", seq.id=sid, order.by=tid, pca_mdl = pca_mdl)
  list(model = mdl, 
       subsets = predict(mdl, dat))
}

#' @export
predict.stagePCACP = function(obj, data) {
  data = as.data.table(data)
  datMat = as.matrix(data[, setdiff(names(data), c(obj$order.by, obj$seq.id)), , with=F])
  data$PC1 = predict(obj$pca_mdl, datMat)[, "PC1"]
  
  cptFn = function(x) {
    mdl = cpt.mean(x, penalty="SIC", method="AMOC", class=FALSE)
    ret = rep("Stage1", length(x))
    if (mdl[["conf.value"]] > .5)
      ret[mdl[["cpt"]]:length(x)] = "Stage2"
    ret
  }
  
  data[, stage := cptFn(.SD[order(get(obj$order.by)), PC1]), by = get(obj$seq.id)]
  list(Stage1 = which(data$stage == "Stage1"),
       Stage2 = which(data$stage == "Stage2"))
}

#' @export
stageThClassifier = function(learner=makeLearner("classif.ksvm"), threshold, predictWindow) {
  assertClass(learner, "RLearnerClassif")
  function(task) {
    td = getTaskDesc(task)
    tid = td$order.by; sid = td$seq.id; tar = td$target
    dat = getTaskData(task)
    # Make classifier task
    clDat = dat[, setdiff(names(dat), c(tid, sid, tar))]
    if (! threshold %between% range(dat[[tar]]))
      stop("Threshold does not divide data!")
    clDat$.Stage = ifelse(dat[[tar]] <= threshold, "Stage2", "Stage1")
    clTask = makeClassifTask(data = clDat, target = ".Stage", positive = "Stage2")
    clMdl = train(learner, clTask)
    
    mdl = makeS3Obj(classes = "stageThClassifier", seq.id=sid, order.by=tid, predictWindow = predictWindow, clLrn = learner, clMdl = clMdl)
    list(model = mdl, 
         subsets = predict(mdl, dat))
  }
}
#' @export
predict.stageThClassifier = function(obj, data) {
  clDat = data[, setdiff(names(data), c(obj$seq.id, obj$order.by))]
  pDat = as.data.table(data[, c(obj$seq.id, obj$order.by)])
  pDat$.stage = predictLearner(.learner = obj$clLrn, .model = obj$clMdl, .newdata = clDat)
  pDat[order(get(obj$order.by)), 
       stage_filtered := {
        tmp = cumsum(.stage == "Stage2")
        ifelse(tmp > obj$predictWindow, "Stage2", "Stage1")
       }, 
       by=get(obj$seq.id)]
  list(Stage1 = which(pDat$stage_filtered == "Stage1"),
       Stage2 = which(pDat$stage_filtered == "Stage2"))
}