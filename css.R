# ============================================================
# css.R - Returns app CSS by reading the installed stylesheet
# ============================================================

#' Read the app CSS from the installed package file
#' @keywords internal
app_css <- function() {
  css_file <- system.file("app/www/styles.css", package = "SeuratExplorer")
  if (!nzchar(css_file)) {
    css_file <- file.path(
      system.file(package = "SeuratExplorer"), "styles.css"
    )
  }
  paste(readLines(css_file, warn = FALSE), collapse = "\n")
}
