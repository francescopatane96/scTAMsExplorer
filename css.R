# ============================================================
# css.R - Returns app CSS by reading the installed stylesheet
# ============================================================

#' Read the app CSS from the installed package file
#' @keywords internal
app_css <- function() {
  # system.file() cerca dentro inst/ del pacchetto installato:
  # se il file è inst/styles.css, passa "www/styles.css"
  # se è inst/styles.css, passa "styles.css"
  css_file <- system.file("/inst/app/www/styles.css", package = "scTAMsExplorer")

  if (!nzchar(css_file)) {
    css_file <- system.file("styles.css", package = "scTAMsExplorer")
  }

  if (!nzchar(css_file) || !file.exists(css_file)) {
    warning("styles.css not found in package; app will render without custom CSS.")
    return("")
  }

  paste(readLines(css_file, warn = FALSE), collapse = "\n")
}
