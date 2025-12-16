# Startup message for scTrends ------------------------------------------------

#' @keywords internal
scTrends_startup_message <- function() {
  v <- utils::packageVersion("scTrends")

  # Unicode-safe symbols (ASCII source, Unicode display)
  ball_filled <- "\u2B22"  # ⬢
  ball_open   <- "\u2B21"  # ⬡

  # Colored symbols line 1
  line_top <- paste0(
    "  ",
    crayon::blue(ball_filled),  "   ",
    crayon::magenta(ball_open), "       ",
    crayon::green(ball_filled), "      ",
    crayon::yellow(ball_open), "\n"
  )

  # ASCII logo: scTrends
  logo <- paste0(
    crayon::silver("   ____      _____                    _      \n"),
    crayon::silver("  / ___| ___|_   _|__ _ __ _ __   ___| |__   \n"),
    crayon::silver("  \\___ \\/ _ \\| |/ _ \\ '__| '_ \\ / __| '_ \\  \n"),
    crayon::silver("   ___) |  __/ |  __/ |  | | | | (__| | | | \n"),
    crayon::silver("  |____/ \\___|_|\\___|_|  |_| |_|\\___|_| |_| \n"),
    crayon::silver("                 scTrends                   \n")
  )

  # Colored symbols line 2
  line_mid <- paste0(
    "  ",
    crayon::green(ball_filled), " ",
    crayon::silver(" single-cell pseudotime gene trends detection "),
    crayon::magenta(ball_open),
    "\n"
  )

  bar <- crayon::silver(paste(rep("-", 60), collapse = ""))

  info <- paste0(
    bar, "\n",
    "scTrends version ", crayon::yellow(v), "\n\n",
    crayon::silver("Bin, tune, classify, and visualize gene expression trends\n"),
    crayon::silver("along single-cell pseudotime.\n\n"),
    crayon::silver("To suppress this message:\n"),
    "  suppressPackageStartupMessages(library(scTrends))\n",
    bar, "\n"
  )

  paste0(line_top, logo, line_mid, "\n", info)
}

.onAttach <- function(libname, pkgname) {

  if (!interactive()) return(invisible(NULL))

  packageStartupMessage(scTrends_startup_message())
}
