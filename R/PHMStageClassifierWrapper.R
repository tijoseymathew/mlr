#' @title PHM with stage classification
#'
#' @description
#' Wraps a PHM learner that initially trains a simple classifier to detect normal/abnormal stage and trains PHM learner only on abnormal data
#' @export
makePHMStageClassifierWrapper = function(learner) {
  checkLearner(learner, "phmregr")
  ps = makeParamSet(
    makeUntypedLearnerParam("stage_learner", tunable = FALSE, default = stagePCAChangePoint),
    makeUntypedLearnerParam("transition_fn", tunable = FALSE, default = stageTransitionHard)
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
                                                  stage_learner = stagePCAChangePoint, # Returns a model and a changepoint data.table with Stage_startT, Stage_endT, and stage label for each seq.id,
                                                  transition_fn = stageTransitionHard, # Takes data table with .newdata seq.id, order.by and output of stage_learner. Returns weights for each instance for each stage
                                                  ...) {
  assertClass(.task, "PHMRegrTask")
  
  task = subsetTask(.task, .subset)
  # Stage detection
  stage_obj = stage_learner(task)
  
  cpDT = stage_obj$cpDT
  td = getTaskDesc(task)
  tid = td$order.by; sid = td$seq.id
  tDT = as.data.table(getTaskData(task, features = c(sid, tid)))
  next.mdls = lapply(unique(cpDT$Stage),
                     function(stg) {
                       ix = tDT[, {
                         t = get(tid)
                         cp = cpDT[Stage==stg & get(sid) == .BY[[1]]]
                         t > cp$Stage_startT & t <= cp$Stage_endT
                       }, by = list(get(sid))]$V1
                       train(learner = .learner$next.learner,
                             task = task,
                             subset = ix,
                             weights = .weights[ix])
                     })
  names(next.mdls) = unique(cpDT$Stage)
  
  m = makeWrappedModel.Learner(learner = .learner, 
                               learner.model = list(next.model = next.mdls,
                                                    stage.model = stage_obj$model,
                                                    seq.id = sid, order.by = tid,
                                                    transition_fn = transition_fn),
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
  sid = mdl$seq.id; tid = mdl$order.by
  
  # Stage detector
  cpDT = predict(mdl$stage.model, .newdata)
  nDT = as.data.table(.newdata[, c(sid, tid)])
  names(nDT) = c("seq.id", "order.by")
  names(cpDT) = c("seq.id", names(cpDT)[-1])
  transMat = mdl$transition_fn(nDT, cpDT)
  
  yMat = matrix(0, nrow=nrow(.newdata), ncol = length(mdl$next.model))
  colnames(yMat) = names(mdl$next.model)
  if (length(intersect(colnames(transMat), colnames(yMat))) != ncol(yMat))
    stop("Something is wrong. All stages not present somewhere")
  for (stg in names(mdl$next.model)) {
    ix = which(transMat[, stg] != 0)
    if (length(ix))
      yMat[ix, stg] = do.call( predictLearner,
                               c( list(.learner = .learner$next.learner,
                                       .model = mdl$next.model[[stg]],
                                       .newdata = .newdata[ix, ]),
                                  args) )$y
  }
  cbind(.newdata[, c(sid, tid)], y=rowSums(yMat*transMat))
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
       cpDT = predict(mdl, dat))
}

#' @export
predict.stagePCACP = function(obj, data) {
  data = as.data.table(data)
  datMat = as.matrix(data[, setdiff(names(data), c(obj$order.by, obj$seq.id)), , with=F])
  data$PC1 = predict(obj$pca_mdl, datMat)[, "PC1"]
  
  cptFn = function(x) {
    mdl = cpt.mean(x, penalty="SIC", method="AMOC", class=FALSE)
    ret = length(x)
    if (!is.na(mdl[["conf.value"]]))
      if (mdl[["conf.value"]] > .5)
        ret = mdl[["cpt"]]
    ret
  }
  
  ret = data[order(get(obj$order.by)), {
    cp = cptFn(PC1)
    t = get(obj$order.by)
    list(c(min(t)+.Machine$double.eps, t[cp]),
         c(t[cp], max(t)+1), #FIXME?
         c("Stage1", "Stage2"))
  }, 
  by = list(get(obj$seq.id))]
  names(ret) = c(obj$seq.id, "Stage_startT", "Stage_endT", "Stage")
  ret
}

#' #' @export
#' stageThClassifier = function(learner=makeLearner("classif.ksvm"), threshold, predictWindow) {
#'   assertClass(learner, "RLearnerClassif")
#'   function(task) {
#'     td = getTaskDesc(task)
#'     tid = td$order.by; sid = td$seq.id; tar = td$target
#'     dat = getTaskData(task)
#'     # Make classifier task
#'     clDat = dat[, setdiff(names(dat), c(tid, sid, tar))]
#'     if (! threshold %between% range(dat[[tar]]))
#'       stop("Threshold does not divide data!")
#'     clDat$.Stage = ifelse(dat[[tar]] <= threshold, "Stage2", "Stage1")
#'     clTask = makeClassifTask(data = clDat, target = ".Stage", positive = "Stage2")
#'     clMdl = train(learner, clTask)
#'     
#'     mdl = makeS3Obj(classes = "stageThClassifier", seq.id=sid, order.by=tid, predictWindow = predictWindow, clLrn = learner, clMdl = clMdl)
#'     list(model = mdl, 
#'          subsets = predict(mdl, dat))
#'   }
#' }
#' #' @export
#' predict.stageThClassifier = function(obj, data) {
#'   clDat = data[, setdiff(names(data), c(obj$seq.id, obj$order.by))]
#'   pDat = as.data.table(data[, c(obj$seq.id, obj$order.by)])
#'   pDat$.stage = predictLearner(.learner = obj$clLrn, .model = obj$clMdl, .newdata = clDat)
#'   pDat[order(get(obj$order.by)), 
#'        stage_filtered := {
#'          tmp = cumsum(.stage == "Stage2")
#'          ifelse(tmp > obj$predictWindow, "Stage2", "Stage1")
#'        }, 
#'        by=get(obj$seq.id)]
#'   list(Stage1 = which(pDat$stage_filtered == "Stage1"),
#'        Stage2 = which(pDat$stage_filtered == "Stage2"))
#' }

#' @export
stageTransitionHard = function(data, cp) {
  ret = lapply(unique(cp$Stage), 
               function(stg) {
                 as.numeric(
                   cp[Stage==stg][data, order.by > Stage_startT & order.by <= Stage_endT, on="seq.id"]
                 )})
  names(ret) = unique(cp$Stage)
  as.matrix(as.data.frame(ret))
}

#' @export
stageTransitionLinear = function(data, cp, width=21L) {
  if (width %% 2 != 1)
    stop("Width should be odd")
  upWts = (1/width) * 1:width
  dwWts = rev(upWts)
  repFn = function(x) {
    d = diff(x)
    for (ix in which(d!=0)) {
      six = max(1, ix-(width %/% 2)); eix = min(length(x), ix+(width %/% 2))
      x[six:eix] = (if (d[ix]==1) upWts else dwWts)[1:(eix-six+1)]
    }
    x 
  }
  linWt = as.data.table(stageTransitionHard(data, cp))
  linWt[, (c("seq.id", "order.by")) := data[, .(seq.id, order.by)]]
  linWt = melt(linWt, id.vars=c("seq.id", "order.by"))
  linWt[order(order.by), value := repFn(value), by=.(seq.id, variable)]
  linWt = dcast(linWt, seq.id+order.by~variable, value.var = "value")
  linWt = linWt[data[, .(seq.id, order.by)], , on=c("seq.id", "order.by")]
  linWt$seq.id = linWt$order.by = NULL
  linWt = as.matrix(linWt)
  sweep(linWt, 1, rowSums(linWt), "/")
}
