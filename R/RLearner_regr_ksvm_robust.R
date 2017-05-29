#' @export
makeRLearner.regr.ksvm.robust = function() {
  makeRLearnerRegr(
    cl = "regr.ksvm.robust",
    package = "kernlab",
    par.set = makeParamSet(
      # FIXME: scaled is not implemented!
      makeLogicalLearnerParam(id = "scaled", default = TRUE),
      makeDiscreteLearnerParam(id = "type", default = "eps-svr", values = c("eps-svr", "nu-svr", "eps-bsvr")),
      makeDiscreteLearnerParam(id = "kernel", default = "rbfdot",
                               values = c("vanilladot", "polydot", "rbfdot", "tanhdot", "laplacedot", "besseldot", "anovadot", "splinedot")),
      makeNumericLearnerParam(id = "C",
                              lower = 0, default = 1, requires = quote(type %in% c("eps-svr", "eps-bsvr"))),
      makeNumericLearnerParam(id = "nu",
                              lower = 0, default = 0.2, requires = quote(type == "nu-svr")),
      makeNumericLearnerParam(id = "epsilon", lower = 0, default = 0.1,
                              requires = quote(type %in% c("eps-svr", "nu-svr", "eps-bsvr"))),
      makeNumericLearnerParam(id = "sigma",
                              lower = 0, requires = quote(kernel %in% c("rbfdot", "anovadot", "besseldot", "laplacedot"))),
      makeIntegerLearnerParam(id = "degree", default = 3L, lower = 1L,
                              requires = quote(kernel %in% c("polydot", "anovadot", "besseldot"))),
      makeNumericLearnerParam(id = "scale", default = 1, lower = 0,
                              requires = quote(kernel %in% c("polydot", "tanhdot"))),
      makeNumericLearnerParam(id = "offset", default = 1,
                              requires = quote(kernel %in% c("polydot", "tanhdot"))),
      makeIntegerLearnerParam(id = "order", default = 1L,
                              requires = quote(kernel == "besseldot")),
      makeNumericLearnerParam(id = "tol", default = 0.001, lower = 0),
      makeLogicalLearnerParam(id = "shrinking", default = TRUE),
      makeLogicalLearnerParam(id = "fit", default = TRUE, tunable = FALSE),
      makeIntegerLearnerParam(id = "cache", default = 40L, lower = 1L),
      makeNumericLearnerParam(id = "gamma", lower = 0, default = 0) # Numerical regularization parameter
    ),
    par.vals = list(fit = FALSE),
    properties = c("numerics", "factors"),
    name = "Support Vector Machines",
    short.name = "ksvm.robust",
    note = "Kernel parameters have to be passed directly and not by using the `kpar` list in `ksvm`. Note that `fit` has been set to `FALSE` by default for speed. Scaled is not implemented.",
    callees = "ksvm.robust"
  )
}

#' @export
trainLearner.regr.ksvm.robust = function(.learner, .task, .subset, .weights = NULL, 
                                         degree, offset, scale, sigma, order, length, lambda, 
                                         scaled = TRUE, gamma = 0, ...) {
  task = subsetTask(.task, .subset)
  scaleMdl = getKSVMScaling(task, scaled)
  task = scaleMdl$task
  
  kFn = getKSVMKFunction(task,
                         degree, offset, scale, sigma, order, length, lambda, 
                         ...)
  
  tr.df = getTaskData(task, target.extra = TRUE)$data
  kMat = kernlab::kernelMatrix(kFn, as.matrix(tr.df))
  kMat = kMat + gamma * diag(nrow = nrow(kMat), ncol = ncol(kMat))
  
  mdl = kernlab::ksvm(kMat, getTaskTargets(task), kernel = "matrix", ...)
  
  svix = kernlab::SVindex(mdl)
  tr.df = tr.df[svix, , drop = FALSE]
  
  list(scale.mdl = scaleMdl$model,
       ksvm.mdl = mdl, 
       kFn = kFn,
       sv.df = tr.df)
}

#' @export
predictLearner.regr.ksvm.robust = function(.learner, .model, .newdata, ...) {
  .newdata = KSVMScaling.x(.model$learner.model$scale.mdl, .newdata)
  kMat = kernlab::kernelMatrix(.model$learner.model$kFn, 
                               as.matrix(.newdata),
                               as.matrix(.model$learner.model$sv.df))
  
  predres = kernlab::predict(.model$learner.model$ksvm.mdl, newdata = kMat, ...)[, 1L]
  KSVMScaling.y(.model$learner.model$scale.mdl, predres)
}


# Additional helper functions ---------------------------------------------

getKSVMScaling = function(task, scaled) {
  tData = getTaskData(task, target.extra = TRUE)
  x = as.matrix(tData$data)
  y = tData$target
  
  x.scale = y.scale = NULL
  if (length(scaled) == 1)
    scaled <- rep(scaled, ncol(x))
  if (any(scaled)) {
    co <- !apply(x[,scaled, drop = FALSE], 2, var)
    if (any(co)) {
      scaled <- rep(FALSE, ncol(x))
      warning(paste("Variable(s)",
                    paste("`",colnames(x[,scaled, drop = FALSE])[co],
                          "'", sep="", collapse=" and "),
                    "constant. Cannot scale data.")
      )
    } else {
      xtmp <- scale(x[,scaled])
      x[,scaled] <- xtmp
      x.scale <- attributes(xtmp)[c("scaled:center","scaled:scale")]
      # if (is.numeric(y)&&(type(ret)!="C-svc"&&type(ret)!="nu-svc"&&type(ret)!="C-bsvc"&&type(ret)!="spoc-svc"&&type(ret)!="kbb-svc")) {
      if (is.numeric(y)) {
        y <- scale(y)
        y.scale <- attributes(y)[c("scaled:center","scaled:scale")]
        y <- as.vector(y)
      }
    }
  }
  tData = data.frame(x)
  tData[[getTaskTargetNames(task)]] = y
  list(task = changeData(task, tData),
       model = list(scaled=scaled, x.scale = x.scale, y.scale = y.scale))
}

KSVMScaling.x = function(model, newdata) {
  if (any(model$scaled))
    newdata[,model$scaled] <-
      scale(newdata[,model$scaled, drop = FALSE],
            center = model$x.scale$"scaled:center", scale = model$x.scale$"scaled:scale")
  newdata
}

KSVMScaling.y = function(model, predres) {
  # if (!is.null(model$y.scale) & !is(newdata,"kernelMatrix") & !is(newdata,"list"))
  if (any(model$scaled))
    if (!is.null(model$y.scale))
      ## return raw values, possibly scaled back
      return(predres * model$y.scale$"scaled:scale" + model$y.scale$"scaled:center")
  return(predres)
}

getKSVMKFunction = function(task, degree, offset, scale, sigma, order, length, lambda, ...) {
  args = list(...)
  kernel = ifelse("kernel" %in% names(args), args$kernel, "rbfdot")
  if("scaled" %in% names(args)) {
    if(args$scaled)
      stop("Not sure how to handle scaling :-(")
  }
  
  td = getTaskDesc(task)
  # Extract kernel function
  kpar = learnerArgsToControl(list, degree, offset, scale, sigma, order, length, lambda)
  if (is.character(kernel)) {
    kernel = match.arg(kernel, c("rbfdot", "polydot", "tanhdot", "vanilladot", "laplacedot", "besseldot", "anovadot", "splinedot", "matrix"))
    
    matchDefault <- function(p, default) ifelse(p %in% names(kpar), kpar[[p]], default)
    
    sigma <- matchDefault("sigma",
                          if (kernel %in% c("rbfdot", "laplacedot") | identical(kpar, "automatic")) {
                            datMat = as.matrix(BBmisc::dropNamed(getTaskData(task),
                                                                 c(td$seq.id,
                                                                   td$order.by,
                                                                   td$target)))
                            mean(kernlab::sigest(datMat, scaled=FALSE)[c(1L, 3L)])
                          } else
                            0.1)
    degree <- matchDefault("degree", 1)
    offset <- matchDefault("offset", 1)
    scale <- matchDefault("scale", 1)
    order <- matchDefault("order", 1)
    kFn = switch (kernel,
                  rbfdot = kernlab::rbfdot(sigma = sigma),
                  polydot = kernlab::polydot(degree = degree, scale = scale, offset = offset),
                  tanhdot = kernlab::tanhdot(scale = scale, offset = offset),
                  vanilladot = kernlab::vanilladot(),
                  laplacedot = kernlab::laplacedot(sigma = sigma),
                  besseldot = kernlab::besseldot(sigma = sigma, order = order, degree = degree),
                  anovadot = kernlab::anovadot(sigma = sigma, degree = degree),
                  splinedot = kernlab::splinedot(),
                  matrix = stopf("kMat interface requires some work :-("), 
                  stopf("Provided kernel %s is not implemented", kernel)
    )
  } else if (is.function(kernel)) 
    kFn = kernel
  kFn
}