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
  ps = makeParamSet(makeFunctionLearnerParam("clust_fn"),
                    makeUntypedLearnerParam("ens_lrn"))
  
  lrn = makeBaseWrapper(id = id,
                        type = "phmregr", next.learner = learner,
                        par.set = ps,
                        learner.subclass = "PHMEnsembleWrapper",
                        model.subclass = "PHMEnsembleModel")
  lrn$properties = learner$properties
  return(lrn)
}

#' @export
print.PHMEnsembleModel = function(x, ...) {
  s = capture.output(print.WrappedModel(x))
  u = sprintf("PHM Ensemble Learner: %s", class(x$learner$next.learner)[1L])
  s = append(s, u, 1L)
  lapply(s, catf)
}

#' @export
trainLearner.PHMEnsembleWrapper = function(.learner, .task, .subset, .weights = NULL, clust_fn, ens_lrn = NULL, ...) {
  
  assertClass(.task, "PHMRegrTask")
  .task = subsetTask(.task, subset = .subset)
  
  args = list(task = .task, learner = .learner, weights = .weights)
  parallelLibrary("mlr", master = FALSE, level = "mlr.ensemble", show.info = FALSE)
  exportMlrOptions(level = "mlr.ensemble")
  clust_model = clust_fn(.task)
  models = parallelMap(doPHMEnsembleTrainIteration, clust_model$clusters, more.args = args, level = "mlr.ensemble")
  ens_mdl = NULL
  tid = .task$task.desc$order.by; sid = .task$task.desc$seq.id
  if (!is.null(ens_lrn)) {
    checkLearner(ens_lrn, type = c("phmregr"))
    p = asMatrixCols(lapply(models, function(m) {
      predict(m, .task, ...)$data$response
    }))
    colnames(p) = paste0("p", 1:ncol(p))
    ens_tsk = makePHMRegrTask(data=cbind(p, getTaskData(.task, features = c(tid, sid))),
                              target = getTaskTargetNames(.task), seq.id = sid, order.by = tid)
    ens_mdl = train(ens_lrn, ens_tsk)
  }
  
  m = makeWrappedModel.Learner(learner = .learner, 
                               learner.model = list(clust_model = clust_model$model,
                                                    ensemble_models = models,
                                                    ens_lrn = ens_lrn,
                                                    ens_mdl = ens_mdl,
                                                    order.by = tid,
                                                    seq.id = sid),
                               task.desc = getTaskDesc(.task),
                               subset = .subset,
                               features = getTaskFeatureNames(.task),
                               factor.levels = mlr:::getTaskFactorLevels(.task),
                               time = 0)
  makeChainModel(next.model = m, cl = "PHMEnsembleWrapperModel")
}

doPHMEnsembleTrainIteration = function(seq.id, task, learner, weights) {
  setSlaveOptions()
  sids = getTaskData(task, features = getTaskDesc(task)$seq.id, target.extra = TRUE)$data[, 1]
  ix = (sids %in% seq.id)
  train(learner$next.learner, task, subset = ix, weights = weights[ix])
}

#' @export
predictLearner.PHMEnsembleWrapper = function(.learner, .model, .newdata, ...) {
  mdl = .model$learner.model$next.model$learner.model
  models = mdl$ensemble_models
  p = asMatrixCols(lapply(models, function(m) {
    predict(m, newdata = .newdata, ...)$data$response
  }))
  if (!is.null(mdl$ens_mdl)) {
    colnames(p) = paste0("p", 1:ncol(p))
    p = cbind(p, .newdata[, c(mdl$seq.id, mdl$order.by)])
    return(predictLearner(mdl$ens_lrn, mdl$ens_mdl, p))
  }
  weights = predict(mdl$clust_model, .newdata)
  ret = predict(models[[1]], newdata = .newdata, ...)$data
  ret$response = NULL
  ret$y = rowSums(p * weights, na.rm = TRUE)
  ret
}

#' Helper function to get predictions from each ensemble members
#' @export
predictPHMEnsemble = function(model, newdata, ...) {
  assertClass(model, "PHMEnsembleModel")
  models = getLearnerModel(model, more.unwrap = FALSE)
  
  ret = predict(models[[1]], newdata, ...)$data
  r_col = which(grepl("^response$", colnames(ret)))
  colnames(ret)[r_col] = "response_1"
  i = 2
  for (mod in models[-1]){
    p = predict(mod, newdata, ...)$data
    colnames(p)[r_col] = paste0("response_", i)
    i = i + 1
    ret = merge(ret, p)
  }
  ret
}

#' Helper function to get prediction from each ensemble learner in a benchmark test
#' WARNING: May take a lot of time, consider passing a small task
#' @export
getBMRPHMEnsemblePredictions = function(bmr, learner.ids = getBMRLearnerIds(bmr), task) {
  assertClass(bmr, "BenchmarkResult")
  assertSubset(learner.ids, getBMRLearnerIds(bmr))
  ret = lapply(learner.ids, function(x) {
    m = bmr$results$cmapss[[x]]$models[[1]]$learner.model$next.model
    p = as.data.table(predictPHMEnsemble(m, task))
    p[, learner.id := x]
  })
  as.data.frame(Reduce(rbind, ret))
}

#' Plot BMRPHMEnsemble predictions as a distribution
#' Only last prediction for each seq.id will be plotted
#' WARNING: May take a lot of time, consider passing a small task
#' @export
plotBMRPHMEnsembleDistribution = function(bmr, learner.ids, task, seq.ids = 1) {
  checkTask(task, "PHMRegrTask")
  td = getTaskDesc(task)
  sid = td$seq.id; tid = td$order.by
  
  t_s.ids = task$env$data[[sid]]
  if (is.numeric(seq.ids)) {
    assertNumeric(seq.ids, lower = 0, upper = 1)
    u_sids = unique(t_s.ids)
    seq.ids = sample(u_sids, round(length(u_sids) * seq.ids))
  }
  task = subsetTask(task, t_s.ids %in% seq.ids)
  plotDT = as.data.table( getBMRPHMEnsemblePredictions(bmr, learner.ids, task) )
  plotDT = plotDT[, .SD[get(tid) == max(get(tid))], by = sid]
  rCols = c("learner.id", "truth", tid, sid, grep("response_[0-9]*", colnames(plotDT), value = TRUE))
  plotDT[[setdiff(colnames(plotDT), rCols)]] = NULL
  plotDT = melt(plotDT, id.var = c("learner.id", tid, sid, "truth"))
  ggplot(plotDT, aes(col = learner.id, group=learner.id))+
    geom_density(aes(value))+
    geom_vline(data = plotDT[, .SD[1], by = sid], aes(xintercept = truth))+
    facet_grid(as.formula(paste0(sid, "~.")))
}

#' Plots the prediction density for one seq.id on BMRPHMEnsemble models
#' @export
plotBMRPHMEnsemblePrediction = function(bmr, learner.ids, task, seq.id, n.regions = 4) {
  checkTask(task, "PHMRegrTask")
  td = getTaskDesc(task)
  sid = td$seq.id; tid = td$order.by
  
  t_s.ids = task$env$data[[sid]]
  if (missing(seq.id))
    seq.id = sample(t_s.ids, 1)
  assertSubset(seq.id, t_s.ids)
  task = subsetTask(task, t_s.ids == seq.id)
  plotDT = as.data.table( getBMRPHMEnsemblePredictions(bmr, learner.ids, task) )
  
  lb_q = seq(0, 0.5, length.out = n.regions+1)[-(n.regions+1)]
  ub_q = seq(0.5, 1, length.out = n.regions+1)[-1]
  rCols = grep("response_[0-9]*", colnames(plotDT), value=TRUE)
  plotDT[, response_mean := apply(.SD, 1, mean), .SDcols = rCols]
  lb_quants = t(apply(plotDT[, rCols, with=FALSE], 1, quantile, probs = lb_q))
  colnames(lb_quants) = 1:n.regions
  ub_quants = t(apply(plotDT[, rCols, with=FALSE], 1, quantile, probs = ub_q))
  colnames(ub_quants) = n.regions:1
  lbDT = melt(cbind(plotDT[, c("learner.id", "time")], lb_quants),
              id.vars = c("learner.id", "time"), value.name = "lb")
  ubDT = melt(cbind(plotDT[, c("learner.id", "time")], ub_quants),
              id.vars = c("learner.id", "time"), value.name = "ub")
  quantDT = merge(lbDT, ubDT)
  
  ggplot(melt(plotDT, measure.vars = c("truth", "response_mean")), aes(x=time))+
    geom_ribbon(data = quantDT, aes(ymin=lb, ymax=ub, alpha = variable), fill="green")+
    scale_alpha_manual(values = seq(0.25, 0.45, length.out = n.regions))+
    geom_line(aes(y=value, col = variable))+
    scale_color_manual(values = c("red", "blue"))+
    guides(alpha = "none")+
    facet_grid(learner.id~.)+
    ggtitle(seq.id)
}

#' Compute accuracy table of bmr PHMEnsemble prediction being within given probability bounds
#' Assumes a Gaussian distribution of predictions
#' @export
getBMRPHMEnsembleAccuracy = function(bmr, learner.ids, task, p_range = 0.9, only_last = TRUE) {
  checkTask(task, "PHMRegrTask")
  assertNumeric(p_range, 0, 1)
  td = getTaskDesc(task)
  sid = td$seq.id; tid = td$order.by
  
  predDT = as.data.table( getBMRPHMEnsemblePredictions(bmr, learner.ids, task) )
  if (only_last)
    predDT = predDT[, .SD[get(tid) == max(get(tid))], by = sid]
  resCols = grep("response_[0-9]*", colnames(predDT), value = TRUE)
  predDT[, c("lb", "ub", "p_val"):= ({
    mat = .SD[, resCols, with = FALSE]
    m = apply(mat, 1, mean)
    s = apply(mat, 1, sd)
    list(qnorm((1-p_range)/2, mean = m, sd = s),
         qnorm(p_range + (1-p_range)/2, mean = m, sd = s),
         dnorm(truth, mean = m, sd = s))
  })]
  predDT[, in_range := (truth >= lb & truth <= ub)]
  as.data.frame( predDT[, list(accuracy=sum(in_range)/.N,
                               llikelihood=sum(log(p_val))), by=learner.id] )
}
