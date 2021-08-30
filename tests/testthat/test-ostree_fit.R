
devtools::load_all()

f1 = ostree_fit(flchain_x, flchain_y, leaf_min_event = 20)



test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
