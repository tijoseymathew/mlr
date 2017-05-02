#' @export
#' @rdname Task
makePHMRegrTask = function(id = deparse(substitute(data)), data, target, seq.id, order.by, weights = NULL, blocking = NULL, fixup.data = "warn", check.data = TRUE) {
  assertString(id)
  assertDataFrame(data)
  assertString(target)
  assertString(seq.id)
  assertString(order.by)
  assertChoice(fixup.data, choices = c("no", "quiet", "warn"))
  assertFlag(check.data)
  
  if (fixup.data != "no") {
    if (is.integer(data[[target]]))
      data[[target]] = as.double(data[[target]])
    if (!is.character(data[[seq.id]]))
      data[[seq.id]] = as.character(data[[seq.id]])
  }
  
  if (check.data) {
    reqCols = c(target, seq.id, order.by)
    w = which.first(reqCols %nin% colnames(data))
    if (length(w) > 0L)
      stopf("Column names of data doesn't contain var: %s", reqCols[w])
    checkTaskData( data, cols = setdiff(colnames(data), reqCols) )
  }
  
  task = makeSupervisedTask("phmregr", data, target, weights, blocking, fixup.data = fixup.data, check.data = FALSE)
  
  task$task.desc = makePHMRegrTaskDesc(id, data, target, seq.id, order.by, weights, blocking)
  addClasses(task, "PHMRegrTask")
}

makePHMRegrTaskDesc = function(id, data, target, seq.id, order.by, weights, blocking) 
{
  td = makeTaskDescInternal("phmregr", id, dropNamed(data, c(seq.id, order.by)), target, weights, blocking)
  td$seq.id = seq.id
  td$n.seqs = length(unique(data[[seq.id]]))
  td$order.by = order.by
  addClasses(td, c("PHMRegrTaskDesc", "SupervisedTaskDesc"))
}

#' @export
print.PHMRegrTask = function(x, print.weights = TRUE, ...) {
  print.Task(x, print.weights)
  td = x$task.desc
  catf("Sequence Id Column: %s", td$seq.id)
  catf("Number of sequences: %i", td$n.seqs)
  catf("Time Column: %s", td$order.by)
}
