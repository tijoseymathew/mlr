#' @title PHM with Health Index Transformation
#'
#' @description
#' Wraps a PHM learner that trains to target a given health index that is scaled back during prediction
#'
#' @template arg_learner
#' @param hi_function [\code{function(task)}]\cr
#'   Function that returns hi: computed HI for task and hi2rulFn: function to compute RUL from HI
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
    makeFunctionLearnerParam("hi_function", tunable = FALSE, default = makeHIFunction())
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
                                     hi_function = makeHIFunction(), ...) {
  assertClass(.task, "PHMRegrTask")
  assertFunction(hi_function, args = "task")
  
  task = subsetTask(.task, .subset)
  hiMdl = hi_function(task)
  assertNumeric(hiMdl$hi, len = getTaskSize(task), any.missing = FALSE)
  assertFunction(hiMdl$hi2rulFn, args = "data")
  
  data = getTaskData(task, target.extra = TRUE)$data
  data$.hi = hiMdl$hi
  td = getTaskDesc(task)
  task = makePHMRegrTask(id = paste0(getTaskId(.task), ".hi"),
                         data = data, target = ".hi",
                         seq.id = td$seq.id, order.by = td$order.by)
  next.mdl = train(.learner$next.learner, task, weights = .weights)
  
  m = makeWrappedModel.Learner(learner = .learner, 
                               learner.model = list(next.model = next.mdl, 
                                                    hi2rulFn = hiMdl$hi2rulFn),
                               task.desc = getTaskDesc(.task),
                               subset = .subset,
                               features = getTaskFeatureNames(.task),
                               factor.levels = mlr:::getTaskFactorLevels(.task),
                               time = 0)
  ret = makeChainModel(next.model = m, cl = "PHMHIModel")
}

#' @export
predictLearner.PHMHIWrapper = function(.learner, .model, .newdata, 
                                       ret.HI = FALSE, # Convinence parameter used for plotting
                                       ...) {
  args = removeFromDots(names(.learner$par.vals), ...)
  
  mdl = .model$learner.model$next.model$learner.model
  ret = do.call( predictLearner, 
                 c( list(.learner = .learner$next.learner, 
                         .model = mdl$next.model,
                         .newdata = .newdata), 
                    args) )
  if (ret.HI)
    return(ret)
  y = mdl$hi2rulFn(ret)
  assertNumeric(y, len = nrow(.newdata), any.missing = FALSE)
  ret$y = y 
  ret
}

# HI Model Plot -----------------------------------------------------------

#' @title PHM HI model exploratory plots
#' 
#' @description 
#' Plots the created synthtic HI and prediction made by learner on a given task
#' 
#' @template arg_learner
#' @template arg_task
#' @param seq.ids
#'   Sequence id to be plotted. Default "all"
#' @export
plotPHMHI = function(learner, task, seq.ids = "all") {
  assertClass(learner, "PHMHIWrapper")
  assertCharacter(seq.ids)
  
  mdl = train(learner, task)
  
  hi_function = ifelse("hi_function" %in% names(getHyperPars(learner)),
                       getHyperPars(learner)$hi_function, makeHIFunction())
  
  pre = predictLearner.PHMHIWrapper(learner, mdl, getTaskData(task, target.extra = TRUE)$data,
                                    ret.HI = TRUE)
  names(pre)[names(pre) == "y"] = "response"
  pre$truth = hi_function(task)$hi
  setDT(pre)
  
  td = getTaskDesc(task)
  tid = td$order.by; sid = td$seq.id
  if (identical(seq.ids, "all"))
    seq.ids = unique(pre[[sid]])
  pre = pre[get(sid) %in% seq.ids]
  pre = melt(pre, id.vars = c(sid, tid), measure.vars = c("truth", "response"),
             variable.name = "type", value.name = "HI")
  
  ggplot(pre, aes_string(x=tid, y="HI", col="type"))+
    geom_line()+
    facet_grid(paste0(sid, "~."))
}

# HI Functions ------------------------------------------------------------

#' @title HI model functions for PHMHIWrapper
#'
#' @description
#' HI from features on a window at start and end of each seq.id
#' 
#' @param type [\code{character(1)}]\cr
#'   rul_scaling: Scales RUL to [0, 1] for each unit
#'   rul_limiting: Limits RUL to multiplier and then scales to [0, 1]
#'   linear: Extrapolates linear model to features fitted on ends for each unit
#'   exponential: Fits a[exp(-b*RUL+c)-exp(c)] on result of linear
#' @param start.n [\code{numeric(1)}]\cr
#'   Number of points at begining to be modelled as HI 1. Default 1
#' @param end.n [\code{numeric(1)}]\cr
#'   Number of points at end to be modelled as HI 0. Default 1
#' @param scaled [\code{logical(1)}]\cr
#'   Should the features be scaled. Default True
#' @param multiplier [\code{numeric(1)}]
#'   Multiplier to be used during prediction. If NULL mean of train RUL will be used. Defaults to NULL
#' @param exp_tol [\code{numeric(1)}]\cr
#'   Warning will be produced if sse error while fitting exponential model is larger than this value
#' @param exp_max_iter [\code{numeric(1)}]\cr
#'   Maximum number of iteration for exponential function with random initializations
#' @export
makeHIFunction = function(type = "rul_scaling", start.n = 1, end.n = 1, scaled = TRUE, multiplier = NULL,
                          exp_tol = 0.3, exp_max_iter = 100) {
  # Per unit HI functions (for one seq.id, in sorted order.by)
  # Both data tables contain RUL column which may or may not be used
  #> trDT : Data table of ends with column hi to be used for training and rest as features
  #> teDT : Data table with features
  rul_scaling = function(trDT, teDT) {
    teDT$RUL/max(teDT$RU)
  }
  rul_limiting = function(trDT, teDT) {
    sapply(teDT$RUL, min, multiplier)/multiplier
  }
  linear = function(trDT, teDT) {
    mdl = lm(hi~.-RUL, trDT)
    predict(mdl, teDT)
  }
  exponential = function(trDT, teDT) {
    x = teDT$RUL
    y = linear(trDT, teDT)
    # Exponential function defintions
    hi = expression(a*(exp(-b*x+c)-exp(c)))
    er = function(p, x, y) {
      a = p$a; b = p$b; c = p$c
      y-eval(hi)
    }
    ja = function(p, x, y) {
      a = p$a; b = p$b; c = p$c
      -c(eval(D(hi, "a")), eval(D(hi, "b")), eval(D(hi, "c")))
    }
    b_sse = Inf
    for (iter in seq(exp_max_iter)) { # Better initialization can avoid this?
      nlsFit = nls.lm(par = list(a=runif(1), b=runif(1),c=runif(1)),
                      fn = er, jac = ja, x = x, y = y,
                      control = nls.lm.control(maxfev = exp_max_iter, maxiter = exp_max_iter))
      sse = sqrt(mean(er(nlsFit$par, x, y)^2))
      if (!is.na(sse) & sse < b_sse) {
        b_mdl = nlsFit
        b_sse = sse
      }
    }
    if (b_sse > exp_tol)
      warning(sprintf("Exponential HI did not converge. sse=%.3f", b_sse))
    a = b_mdl$par$a; b = b_mdl$par$b; c = b_mdl$par$c
    eval(hi)
  }
  hiFns = list(rul_scaling = rul_scaling, rul_limiting = rul_limiting,
               linear = linear, exponential = exponential)
  assertChoice(type, names(hiFns))
  
  # Returned function
  function(task) {
    td = getTaskDesc(task)
    tid = td$order.by; sid = td$seq.id; tar = td$target
    
    if (is.null(multiplier))
      multiplier = mean(getTaskTargets(task))
    
    dat = as.data.table(getTaskData(task))
    names(dat)[names(dat) == tar] = "RUL"
    if (scaled) {
      cols = setdiff(names(dat), c(tid, sid, "RUL"))
      dat[, (cols) := lapply(.SD, scale), .SDcols = cols]
    }
    
    dat[order(get(tid)), hi := {
      # FIXME: Clip start.n and end.n
      marginDT = rbind(head(.SD, start.n), tail(.SD, end.n))
      marginDT$hi = rep(c(1, 0), times = c(start.n, end.n))
      marginDT[[tid]] = NULL
      hiFns[[type]](marginDT, .SD)
    }, by = sid]
    list(hi = dat$hi, hi2rulFn = function(data) multiplier * data$y)
  }
}