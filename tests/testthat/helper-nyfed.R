# The New York Fed's published ACM series, for the tests that compare against
# it.
#
# It is downloaded, never bundled: its terms are the New York Fed's, not the
# package's MIT licence. So these tests run only off CRAN, and skip rather
# than fail when the download does not work. The download is made once per
# session and shared. The sample ends where gsw_monthly does, so the two
# share a date index.
nyfed_acm <- local({
  cached <- NULL
  function() {
    skip_on_cran()
    skip_if_not_installed("readxl")
    if (is.null(cached)) {
      cached <<- tryCatch(
        acm_published(end = "2024-12-31"),
        error = function(e) conditionMessage(e)
      )
    }
    if (is.character(cached)) {
      skip(paste("New York Fed ACM series unavailable:", cached))
    }
    cached
  }
})
