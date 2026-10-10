test_that("gini0 works", {
  income = c(10, 20, 30, 40, 100)
  r = gini0(income)
  expect_true(r >= 0 & r <= 1)
  expect_equal(gini0(c(10,10,10,10)), 0, tolerance=1e-8)
})

test_that("gini works", {
  y = c(10, 20, 30, 40, 100)
  pop = c(100, 150, 200, 250, 300)
  r = gini(y, pop)
  expect_true(r >= 0 & r <= 1)
})

test_that("theil0 works", {
  r = theil0(c(10, 10, 10))
  expect_equal(r, 0, tolerance=1e-8)
})

test_that("theil works", {
  y = c(10, 20, 30)
  pop = c(5, 5, 5)
  r = theil(y, pop)
  expect_true(r >= 0)
})

test_that("theil0_g works", {
  data = data.frame(g=c("A","A","B","B","B"), y=c(10,10,8,6,6))
  r = theil0_g(data, "g", "y")
  expect_type(r, "list")
  expect_named(r, c("theil", "ratio"))
})

test_that("theil_g works", {
  data = data.frame(g=c("A","A","B"), y=c(10,20,30), pop=c(3,2,5))
  r = theil_g(data, "g", "y", "pop")
  expect_type(r, "list")
})

test_that("inequality input validation", {
  expect_error(gini0("a"), "must be a numeric vector")
  expect_error(gini0(c(-1, 2)), "non-negative")
  expect_error(gini(c(1,2), c(3,4,5)), "same length")
  expect_error(theil_g(1:5, "g", "y", "p"), "must be a data frame")
})

# Regression: theil_g2_cross must decompose the ACTUAL total Theil exactly
# (within-group2 weights are income shares, not population shares).
test_that("theil_g2_cross decomposes the true total Theil", {
  d = data.frame(
    industry = c("A","A","A","A","B","B","B","B"),
    area = c("East","East","West","West","East","East","West","West"),
    province = c("Shanghai","Beijing","Sichuan","Yunnan",
                 "Shanghai","Beijing","Sichuan","Yunnan"),
    avg_wage = c(500, 400, 80, 60, 300, 250, 50, 40),
    emp_num = c(100, 80, 90, 70, 120, 100, 80, 60)
  )
  r = theil_g2_cross(d, "industry", "area", "avg_wage", "emp_num")
  # T equals the direct Theil of the raw rows
  expect_equal(unname(r$theil["T"]), theil(d$avg_wage, d$emp_num),
               tolerance = 1e-10)
  # and the three components sum to T exactly
  expect_equal(unname(r$theil["T"]),
               unname(r$theil["Tb"] + r$theil["Tw_b"] + r$theil["Tw_w"]),
               tolerance = 1e-10)
})

test_that("theil_g2_nest decomposes the true total Theil", {
  d = data.frame(
    province = c("A","A","A","A","B","B"),
    city = c("A1","A1","A2","A2","B1","B1"),
    industry = c("Manu","Serv","Manu","Serv","Manu","Serv"),
    y = c(50000, 45000, 60000, 55000, 70000, 65000),
    pop = c(10000, 8000, 15000, 12000, 10000, 8000)
  )
  r = theil_g2_nest(d, "province", "city", "y", "pop")
  expect_equal(unname(r$theil["T"]), theil(d$y, d$pop), tolerance = 1e-10)
  expect_equal(unname(r$theil["T"]),
               unname(r$theil["Tb_g1"] + r$theil["Tb_g2"] + r$theil["Tw"]),
               tolerance = 1e-10)
})

test_that("theil_g2_nest: ratio names are unique and its blocks add up", {
  d = data.frame(
    province = c("A","A","A","A","B","B"),
    city = c("A1","A1","A2","A2","B1","B1"),
    industry = c("Manu","Serv","Manu","Serv","Manu","Serv"),
    y = c(50000, 45000, 60000, 55000, 70000, 65000),
    pop = c(10000, 8000, 15000, 12000, 10000, 8000)
  )
  r = theil_g2_nest(d, "province", "city", "y", "pop")

  # The per-group1 total-within block used to reuse the group1 names, which
  # made it unreachable by name; it is now suffixed ".within".
  expect_false(anyDuplicated(names(r$ratio)) > 0)
  expect_true(any(grepl("\\.within$", names(r$ratio))))

  # theil and ratio differ in length by design (see @return); the three
  # top-level ratios must still sum to 1
  expect_equal(length(r$ratio) - length(r$theil), 2L)
  expect_equal(unname(r$ratio["Tb_g1"] + r$ratio["Tb_g2"] + r$ratio["Tw"]),
               1, tolerance = 1e-10)

  # the per-group1 total-within block is the group1 collapse of the
  # per-group1/group2 within block -- per group, not the global total
  within_total = sum(r$ratio[grepl("\\.within$", names(r$ratio))])
  expect_equal(within_total, unname(r$ratio["Tw"]), tolerance = 1e-10)

  cells = r$ratio[grepl("^[AB]_", names(r$ratio))]
  expect_equal(unname(r$ratio["A.within"]),
               sum(cells[grepl("^A_", names(cells))]), tolerance = 1e-10)
  expect_equal(unname(r$ratio["B.within"]),
               sum(cells[grepl("^B_", names(cells))]), tolerance = 1e-10)
})
