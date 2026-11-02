#--------------------------------------------------------------------------------------
#
#    Render an Rmd report with params given as key=value arguments
#
#--------------------------------------------------------------------------------------
#
# Usage:     Rscript scripts/repl_sc_render_report.R <rmd> <output_file> key1=value1 key2=value2 ...
#
# Purpose:   Wrapper around rmarkdown::render() for report rules whose params are too
#            long for `Rscript -e` (R caps -e expressions at 10,000 bytes), e.g.
#            comma-joined file lists. Each key must be declared in the Rmd YAML
#            header. Values are split on the first "=" only, so they may contain "=".
#
# Notes:     output_file is relative to the Rmd's directory, as with a plain
#            rmarkdown::render() call.
#
#--------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: repl_sc_render_report.R <rmd> <output_file> [key=value ...]")

rmd         <- args[1]
output_file <- args[2]
kv          <- args[-(1:2)]

if (length(kv) > 0 && !all(grepl("=", kv, fixed = TRUE))) {
  stop("Params must be key=value: ", paste(kv[!grepl("=", kv, fixed = TRUE)], collapse = ", "))
}

render_params <- setNames(as.list(sub("^[^=]*=", "", kv)), sub("=.*$", "", kv))

message("Rendering ", rmd, " -> ", output_file)
message("Params: ", paste(names(render_params), collapse = ", "))

rmarkdown::render(rmd, output_file = output_file, params = render_params)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
