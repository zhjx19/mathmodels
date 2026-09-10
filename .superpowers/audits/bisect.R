lines = readLines(".superpowers/audits/probe-pred_grey.R")
env = new.env(parent = globalenv())
for (i in seq_along(lines)) {
  frag = paste(lines[1:i], collapse = "\n")
  ok = tryCatch({ parse(text = frag); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) next   # incomplete fragment (mid-block), skip
  res = tryCatch({ eval(parse(text = frag), envir = env); "ok" },
                 error = function(e) paste0("ERR@", i, ": ", conditionMessage(e)))
  if (res != "ok") {
    cat(res, "\n")
    cat("LINE", i, ":", lines[i], "\n")
    break
  }
}
