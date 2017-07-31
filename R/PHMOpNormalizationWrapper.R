#' @title PHM operational mode normalization
#'
#' @description
#' Normalizes features with respect to different operational modes. 
#' Target is also optionally scaled and re-scaled while prediction (not w.r.t. operational modes)
#' op.feature is removed from task for training with next learner
#'
#' @template arg_learner
#' @param op.feature [\code{character(1)}]\cr
#'   FIXME: Add operational feature to PHM task, add option in learners to incorporate them or not
#'   Feature that determines the operating mode
#' @param method [\code{character(1)}]\cr
#'   Normalizing method. Available are:\cr
#'   \dQuote{center}: Subtract mean.\cr
#'   \dQuote{scale}: Divide by standard deviation.\cr
#'   \dQuote{standardize}: Center and scale.\cr
#'   \dQuote{range}: Scale to a given range.\cr
#' @param range [\code{numeric(2)}]\cr
#'   Range for method \dQuote{range}.
#'   Default is \code{c(0,1)}.
#' @param on.constant [\code{character(1)}]\cr
#'   How should constant vectors be treated? Only used, of \dQuote{method != center},
#'   since this methods does not fail for constant vectors. Possible actions are:\cr
#'   \dQuote{quiet}: Depending on the method, treat them quietly:\cr
#'     \dQuote{scale}: No division by standard deviation is done, input values.
#'        will be returned untouched.\cr
#'     \dQuote{standardize}: Only the mean is subtracted, no division is done.\cr
#'     \dQuote{range}: All values are mapped to the mean of the given range.\cr
#'   \dQuote{warn}: Same behaviour as \dQuote{quiet}, but print a warning message.\cr
#'   \dQuote{stop}: Stop with an error.\cr
#' @param normalize.target [\code{logical(1)}]\cr
#'   Should the target be also normalized
#' @template ret_learner
#' @export
#' @family wrapper
makePHMOpNormalizationWrapper = function(learner) {
  checkLearner(learner, "phmregr")
  ps = makeParamSet(
    makeUntypedLearnerParam("op.feature"),
    makeDiscreteLearnerParam("method", values = c("center", "scale", "standardize", "range"), default = "standardize"),
    makeNumericVectorLearnerParam("range", len = 2, default = c(0, 1)),
    makeDiscreteLearnerParam("on.constant", values = c("quiet", "warn", "stop"), default = "quiet"),
    makeLogicalLearnerParam("normalize.target", default = TRUE)
  )
  lrn = makeBaseWrapper(id = paste(learner$id, "op_normalized", sep = "."),
                        type = "phmregr", next.learner = learner, 
                        par.set = ps,
                        learner.subclass = "PHMOpNormalizationWrapper",
                        model.subclass = "PHMOpNormalizationModel")
  lrn$properties = unique(learner$properties, "factors")
  return(lrn)
}

#' @export
trainLearner.PHMOpNormalizationWrapper = function(.learner, .task, .subset, .weights = NULL, 
                                                  op.feature, method = "standardize", range = c(0,1), on.constant = "quiet", normalize.target = TRUE,
                                                  ...) {
  assertClass(.task, "PHMRegrTask")
  assertSubset(op.feature, getTaskFeatureNames(.task))
  task = subsetTask(.task, .subset)
  
  td = task$task.desc
  sid = td$seq.id; tid = td$order.by; tar = td$target
  data = getTaskData(task)
  norm_cols = setdiff(getTaskFeatureNames(task), c(sid, tid, op.feature))
  assertFactor(data[[op.feature]])
  norm_models = list()
  for (col in norm_cols) {
    norm_models[[col]] = list()
    for (op.mode in levels(data[[op.feature]])) {
      ix = which(data[[op.feature]] == op.mode)
      norm_models[[col]][[op.mode]] = vectorScaleModel(data[[col]][ix], method = method, range = range, on.constant = on.constant)
      data[[col]][ix] = predict(norm_models[[col]][[op.mode]], data[[col]][ix])
    }
  }
  if (normalize.target) {
    norm_models[[tar]] = vectorScaleModel(data[[tar]], method = method, range = range, on.constant = on.constant)
    data[[tar]] = predict(norm_models[[tar]], data[[tar]])
  }
  data[[op.feature]] = NULL
  task = changeData(task, data = data, weights = .weights[.subset])
  n.mdl = train(.learner$next.learner, task, weights = .weights)
  m = makeWrappedModel.Learner(learner = .learner, 
                               learner.model = list(next.model = n.mdl, 
                                                    norm.models = norm_models,
                                                    op.feature = op.feature,
                                                    target = tar),
                               task.desc = getTaskDesc(.task),
                               subset = .subset,
                               features = getTaskFeatureNames(.task),
                               factor.levels = mlr:::getTaskFactorLevels(.task),
                               time = 0)
  makeChainModel(next.model = m, 
                 cl = "PHMOpNormalizationModel")
}

#' @export
predictLearner.PHMOpNormalizationWrapper = function(.learner, .model, .newdata, ...) {
  args = removeFromDots(names(.learner$par.vals), ...)
  
  mdl = .model$learner.model$next.model$learner.model
  for (col in setdiff(names(mdl$norm.models), mdl$target))
    for (op.mode in levels(.newdata[[mdl$op.feature]])) {
      ix = which(.newdata[[mdl$op.feature]] == op.mode)
      .newdata[[col]][ix] = predict(mdl$norm.models[[col]][[op.mode]], .newdata[[col]][ix])
    }
  .newdata[[mdl$op.feature]] = NULL
  ret = do.call(predictLearner, 
                c( list(.learner = .learner$next.learner, 
                        .model = mdl$next.model, 
                        .newdata = .newdata),
                   args)
  )
  if (mdl$target %in% names(mdl$norm.models)) {
    a = mdl$norm.models[[mdl$target]]$a
    b = mdl$norm.models[[mdl$target]]$b
    c = mdl$norm.models[[mdl$target]]$c
    ret$y = (ret$y - c) * b + a
  }
  ret
}

# Returns parameters a,b, c that can be used to scale x as
# (x-a)/b+c
vectorScaleModel = function(x, method, range, on.constant) {
  # is x a constant vector?
  ret = if (length(unique(x[!is.na(x)])) == 1L) {
    switch(on.constant,
           warn = warning("Constant vector in normalization."),
           stop = stop("Constant vector in normalization."))
    switch(method,
           center = list(a = mean(x, na.rm = TRUE), b=1, c=0),
           range = list(a = unique(x[!is.na(x)]) - mean(range), b=1, c=0),
           standardize = list(a = mean(x, na.rm = TRUE), b=1, c=0),
           scale = list(a = 0, b=1, c=0)
    )
  } else {
    switch(method,
           range = list(a = min(x, na.rm = TRUE), b = diff(range(x, na.rm = TRUE)) / diff(range), c = range[1L]),
           standardize = list(a = mean(x, na.rm = TRUE), b = sd(x, na.rm = TRUE), c = 0),
           center = list(a = mean(x, na.rm = TRUE), b = 1, c = 0),
           scale = list(a = 0, b = sd(x, na.rm = TRUE), c = 0)
    )
  }
  makeS3Obj(classes = "vectorScaleModel", a = ret$a, b = ret$b, c = ret$c)
}

predict.vectorScaleModel = function(obj, x) (x - obj$a)/obj$b + obj$c

# FIXME: Above FIXME should eliminate necessity for this hack
getLearnerProperties.PHMOpNormalizationWrapper = function(learner) {
  props = intersect(listLearnerProperties(learner$type), getLearnerProperties(learner$next.learner))
  union(props, "factors")
}