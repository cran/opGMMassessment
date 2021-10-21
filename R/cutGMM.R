# Cuts a vector into separate groups according to input breaks
#' @importFrom methods is
cutGMM <- function(x, breaks, right = TRUE) {
  sizeX <- function(A) {
    TA <- is(A)
    TypeOfA <- TA[2]
    Result <- c(NaN, NaN)
    if (TypeOfA == "vector") {
      Result <- c(length(A), 1)
    }
    if (TypeOfA == "array") {
      Result <- dim(A)
    }
    return(Result)
  }

  if (sizeX(x)[2] > 1 | length(x) == 0)
    stop("CutGMM: x should be a single vector of length > 0.")
  if (length(breaks) > 0) {
    if (hasArg(right) == FALSE) {
      right <- TRUE
    }
    breaks <- sort(breaks)
    if (right == TRUE) {
      if (length(breaks) > 1) {
        compMat <- vapply(x, function(x) x > breaks, logical(length(breaks)))
        ClassesGMM <- colSums(compMat) + 1
      } else ClassesGMM <- as.numeric(x > breaks) + 1
    } else {
      if (length(breaks) > 1) {
        compMat <- vapply(x, function(x) x >= breaks, logical(length(breaks)))
        ClassesGMM <- colSums(compMat) + 1
      } else ClassesGMM <- as.numeric(x >= breaks) + 1
    }
  } else {
    warning("cutGMM: No breaks provided. Assuming one single class.", call. = FALSE)
  }

  # Return results
  return(ClassesGMM)
}

