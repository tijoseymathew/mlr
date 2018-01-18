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
#' @param hi_functions
#'   Either a named list of functions or "all" for all functions in makeHIFunction
#' @param seq.ids [\code{character}]\cr
#'   Sequence id to be plotted. Default "all"
#' @export
plotPHMHI = function(learner, task, hi_functions = "all", seq.ids = "all") {
  assertClass(learner, "PHMHIWrapper")
  assertCharacter(seq.ids)
  if (identical(hi_functions, "all"))
    hi_functions = list(rul_scaling = makeHIFunction("rul_scaling"),
                        rul_limiting = makeHIFunction("rul_limiting"),
                        linear = makeHIFunction("linear"),
                        exponential = makeHIFunction("exponential"))
  if ("hi_function" %in% names(getHyperPars(learner)))
    hi_functions[["learner.par"]] = getHyperPars(learner)$hi_function
  assertList(hi_functions, any.missing = FALSE, types = "function")
  
  td = getTaskDesc(task)
  tid = td$order.by; sid = td$seq.id  
  pltDT = as.data.table( getTaskData(task, features = c(sid, tid),
                                     target.extra = TRUE)$data )
  
  for (ix in seq_along(hi_functions)) {
    fn = hi_functions[[ix]]; nm = names(hi_functions)[ix]
    pltDT[[paste0(nm, "_computed")]] = fn(task)$hi
    
    learner = setHyperPars(learner, hi_function = fn)
    mdl = train(learner, task)
    pre = predictLearner.PHMHIWrapper(learner, mdl, getTaskData(task, target.extra = TRUE)$data,
                                      ret.HI = TRUE)
    pltDT[[paste0(nm, "_predicted")]] = pre$y
  }
  
  if (identical(seq.ids, "all"))
    seq.ids = unique(pltDT[[sid]])
  pltDT = pltDT[get(sid) %in% seq.ids]
  pltDT = melt(pltDT, id.vars = c(sid, tid), variable.name = "HI_Function", value.name = "HI")
  sp = strsplit(as.character(pltDT$HI_Function), "_")
  pltDT$HI_Function = sapply(sp, function(sp) paste0(sp[-length(sp)], collapse = "_"))
  pltDT$type = sapply(sp, function(sp) sp[length(sp)])
  
  ggplot(pltDT, aes_string(x=tid, y="HI", col="HI_Function"))+
    geom_line(aes(linetype = type))+
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
                          exp_tol = 0.3, rPredict = "multiplier") {
  # Per unit HI functions (for one seq.id, in sorted order.by)
  # Both data tables contain RUL column which may or may not be used
  #> trDT : Data table of ends with column hi to be used for training and rest as features
  #> teDT : Data table with features
  rul_scaling = function(trDT, teDT, multiplier) {
    teDT$RUL/max(teDT$RUL)
  }
  rul_limiting = function(trDT, teDT, multiplier) {
    sapply(teDT$RUL, min, multiplier)/multiplier
  }
  linear = function(trDT, teDT, multiplier) {
    mdl = lm(hi~.-RUL, trDT)
    predict(mdl, teDT)
  }
  exponential = function(trDT, teDT, multiplier) {
    x = teDT$.time
    y = linear(trDT, teDT)
    expFit(x[order(x)], y[order(x)])$value[rank(x)]
  }
  expFit = function(x, y) { # x, y should be ordered for filter!
    hi = expression(k*(exp(c)-exp(b*x-c)))
    er = function(p, x, y) {
      k = p$k; b = p$b; c = p$c
      y-eval(hi)
    }
    ja = function(p, x, y) {
      k = p$k; b = p$b; c = p$c
      -c(eval(D(hi, "k")), eval(D(hi, "b")), eval(D(hi, "c")))
    }
    nlsFit = nls.lm(par = list(k=runif(1, 0, .3), b=1/(0.39*max(x)), c=1),
                    fn = er, jac = ja, x = x, y = y, control = nls.lm.control(maxiter = 1000))
    sse = sqrt(mean(er(nlsFit$par, x, y)^2))
    if (sse > exp_tol) {
      if (length(y) > 7 ) {
        y = filter(y, filter = rep(1/7, 7), method = "convolution", sides = 2)
        y[1:3] = 1
        y[(length(y)-2):length(y)] = 0
      }
      nlsFit = nls.lm(par = list(k=runif(1, 0, .3), b=1/(0.39*max(x)), c=1),
                      fn = er, jac = ja, x = x, y = y, control = nls.lm.control(maxiter = 1000))
      sseN = sqrt(mean(er(nlsFit$par, x, y)^2))
      if (sseN > exp_tol) {
        warning(sprintf("Exponential HI did not converge even after filtering. sse_pre=%.3f, sse_post=%.3f", sse, sseN)) 
      } else
        warning(sprintf("Filtering performed for Exponential HI to converge. sse_pre=%.3f, sse_post=%.3f", sse, sseN))
    }
    k = nlsFit$par$k; b = nlsFit$par$b; c = nlsFit$par$c
    list(value=eval(hi), x0 = 2*c/b, sse = sse)
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
    
    dat[, hi := {
      # FIXME: Clip start.n and end.n
      t = .SD[order(get(tid))]
      marginDT = rbind(head(t, start.n), tail(t, end.n))
      marginDT$hi = rep(c(1, 0), times = c(start.n, end.n))
      marginDT[[tid]] = NULL
      te = copy(.SD)
      names(te)[tid == names(te)] = ".time"
      hiFns[[type]](marginDT, te, multiplier)
    }, by = sid]
    hi2rulFn = function(data) multiplier * data$y
    if (type == "exponential" & rPredict == "estimate") {
      hi2rulFn = function(data) {
        data = as.data.table(data)
        data[order(get(tid)), {
          hiMdl = expFit(get(tid), y)
          if (hiMdl$sse > exp_tol) {
            y * multiplier
          } else
            hiMdl$x0 - get(tid)
        }, by = sid]$V1
      }
    }
    list(hi = dat$hi, hi2rulFn = hi2rulFn)
  }
}
