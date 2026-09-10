# Independent audit probe: eval_preprocess.R
devtools::load_all("E:/Github/mathmodels")
options(warn = -1)
P = 0; F = 0
check = function(name, ok, detail = "") {
  if (isTRUE(ok)) { P <<- P + 1; cat("PASS:", name, "-", detail, "\n") }
  else { F <<- F + 1; cat("FAIL:", name, "-", detail, "\n") }
}
x = c(4, 1, NA, 5, 8)

check("standardize vs base::scale",
      identical(standardize(x), as.vector(base::scale(x))))
check("standardize center/scale toggles",
      identical(standardize(x, scale = FALSE), as.vector(base::scale(x, scale = FALSE))))

check("normalize L2",
      isTRUE(all.equal(normalize(c(3,4)), c(3,4)/5, tolerance = 1e-12)))
check("normalize preserves NA", is.na(normalize(x)[3]))

# rescale
check("rescale + maps to [a,b]",
      isTRUE(all.equal(rescale(x), (x - min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)), tolerance=1e-12)))
check("rescale - reverses",
      isTRUE(all.equal(rescale(x, "-"), 1 - (x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)), tolerance=1e-12)))
check("rescale custom range", {
  r = rescale(c(0,10), a=2, b=6); isTRUE(all.equal(r, c(2,6)))
})

# rescale_middle
rm = rescale_middle(6:9, 7)
check("rescale_middle optimal=1, farthest=0",
      abs(rm[which(6:9==7)] - 1) < 1e-12 && abs(min(rm)) < 1e-12,
      sprintf("vals=(%s)", paste(round(rm,3),collapse=",")))

# rescale_interval
ri = rescale_interval(c(35.2,36.0,36.5,37.0,37.5,38.4), 36, 37)
check("rescale_interval inside=1", all(ri[c(2,3,4)] == 1))
M = max(36-35.2, 38.4-37)
check("rescale_interval outside linear",
      abs(ri[1] - (1 - (36-35.2)/M)) < 1e-12 && abs(ri[6] - (1-(38.4-37)/M)) < 1e-12)

# rescale_extreme: CODE is "+": x/max, "-": min/x (standard GRA)
xe = c(2, 4, 8)
check("rescale_extreme + = x/max", isTRUE(all.equal(rescale_extreme(xe, "+"), c(0.25,0.5,1))))
check("rescale_extreme - = min/x", isTRUE(all.equal(rescale_extreme(xe, "-"), c(1,0.5,0.25))))

# rescale_initial / mean
check("rescale_initial + = x/x1", isTRUE(all.equal(rescale_initial(c(2,4,8)), c(1,2,4))))
check("rescale_initial - = x1/x", isTRUE(all.equal(rescale_initial(c(2,4,8), "-"), c(1,0.5,0.25))))
check("rescale_mean = x/mean", isTRUE(all.equal(rescale_mean(c(1,2,3)), c(1,2,3)/2)))

# to_positive
check("to_positive minmax = max - x", isTRUE(all.equal(to_positive(c(1,3,7)), c(6,4,0))))
check("to_positive reciprocal = 1/x", isTRUE(all.equal(to_positive(c(1,2,4), "reciprocal"), c(1,0.5,0.25))))

# edge behavior
check("standardize zero-variance stops",
      inherits(tryCatch(standardize(rep(5,4)), error=\(e) e), "error"))
check("normalize zero-norm stops",
      inherits(tryCatch(normalize(c(0,0,0)), error=\(e) e), "error"))
check("rescale zero-range stops",
      inherits(tryCatch(rescale(rep(2,3)), error=\(e) e), "error"))

cat(sprintf("\n=== eval_preprocess audit: PASS %d / FAIL %d ===\n", P, F))
