test_that("package metadata is complete", {
  desc <- utils::packageDescription("termpremia")

  expect_identical(desc$Package, "termpremia")
  expect_match(desc$License, "MIT")
  expect_true(nzchar(desc$Title))
  expect_true(nzchar(desc$Description))
})
