# Dependency discovery and optional project-local installation.
# Defining this function performs no installation or other side effect.

cp_ensure_dependencies <- function(cfg) {
  # Keep a working package already visible in the current R session.
  if (requireNamespace("nloptr", quietly = TRUE)) return(invisible(TRUE))

  script_dir <- get0(".cp_script_dir", envir = environment(cp_ensure_dependencies),
                     inherits = TRUE, ifnotfound = NULL)
  if (!is.character(script_dir) || length(script_dir) != 1L ||
      is.na(script_dir) || !nzchar(script_dir)) {
    script_dir <- getwd()
  }
  script_dir <- normalizePath(path.expand(script_dir), winslash = "/", mustWork = FALSE)

  project_lib <- cfg$dependency_lib
  if (!is.character(project_lib) || length(project_lib) != 1L ||
      is.na(project_lib) || !nzchar(project_lib)) {
    version_minor <- strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L]
    version_tag <- paste(R.version$major, version_minor, sep = ".")
    project_lib <- file.path(script_dir, "CP_QMLE_R_library", version_tag,
                             R.version$platform)
  }
  project_lib <- normalizePath(path.expand(project_lib), winslash = "/", mustWork = FALSE)
  repos <- cfg$dependency_repos
  if (!is.character(repos) || !length(repos) || anyNA(repos) || any(!nzchar(repos))) {
    repos <- "https://cloud.r-project.org"
  }

  extra_lib <- Sys.getenv("CP_QMLE_R_LIB", unset = "")
  candidates <- c(project_lib, extra_lib,
                  file.path(script_dir, "R_lib"),
                  file.path(script_dir, "..", "work", "R_lib"),
                  file.path(getwd(), "R_lib"),
                  file.path(getwd(), "work", "R_lib"))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates <- unique(vapply(candidates, function(path) {
    normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
  }, character(1L)))
  existing <- candidates[dir.exists(candidates)]
  original_paths <- .libPaths()
  if (length(existing)) {
    # Try each candidate first. This also avoids an incompatible package in
    # an earlier library masking a usable build in a later project library.
    for (library_path in existing) {
      .libPaths(unique(c(library_path, existing, original_paths)))
      if (requireNamespace("nloptr", quietly = TRUE)) {
        message("Using nloptr from an existing R library: ",
                dirname(find.package("nloptr")))
        return(invisible(TRUE))
      }
    }
  }

  lib_literal <- encodeString(project_lib, quote = '"')
  repos_literal <- paste(deparse(repos), collapse = "")
  install_command <- paste0('install.packages("nloptr", lib = ', lib_literal,
                             ", repos = ", repos_literal, ")")
  library_command <- paste0(".libPaths(c(", lib_literal, ", .libPaths()))")
  instructions <- paste0("\nProject library: ", project_lib,
                          "\nRun in the R console:\n  ", install_command,
                          "\n  ", library_command)

  auto_install <- cfg$auto_install_dependencies
  if (is.null(auto_install)) auto_install <- TRUE
  if (!isTRUE(auto_install)) {
    stop("Package 'nloptr' is not available in this R session.", instructions,
         "\nThen rerun the script, or enable CONFIG$auto_install_dependencies.",
         call. = FALSE)
  }

  # .libPaths() silently ignores nonexistent folders, so create it first.
  if (!dir.exists(project_lib)) {
    dir.create(project_lib, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(project_lib)) {
    stop("Could not create the project R library: ", project_lib,
         "\nSet CONFIG$dependency_lib to a writable project folder.",
         instructions, call. = FALSE)
  }
  .libPaths(unique(c(project_lib, existing, original_paths)))
  message("Installing the missing package 'nloptr' into the project library: ", project_lib)
  install_error <- NULL
  tryCatch(
    utils::install.packages("nloptr", lib = project_lib, repos = repos),
    error = function(e) install_error <<- conditionMessage(e)
  )
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    detail <- if (is.null(install_error)) "" else paste0("\nInstallation error: ", install_error)
    stop("The installation did not make package 'nloptr' available.", detail,
         "\nCheck the CRAN connection and whether this library is writable.",
         instructions, "\nAfter installation succeeds, rerun the script.",
         call. = FALSE)
  }
  invisible(TRUE)
}
