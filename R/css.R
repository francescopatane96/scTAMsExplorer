# ============================================================
# css.R - Returns app CSS by reading the installed stylesheet
# ============================================================

#' Read the app CSS from the installed package file
#' @keywords internal
app_css <- function() {
  css_file <- system.file("styles.css", package = "scTAMsExplorer")
  if (!nzchar(css_file)) {
    css_file <- file.path(
      system.file(package = "scTAMsExplorer"), "app", "www", "styles.css"
    )
  }
  paste(readLines(css_file, warn = FALSE), collapse = "\n")
}
