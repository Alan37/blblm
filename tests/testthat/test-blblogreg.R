test_that("blblogreg() and its affiliated functions works", {
  fit <- blblogreg(Direction ~ Lag1 + Lag2 + Volume, data = ISLR::Smarket, m = 3, B = 100)
  expect_s3_class(fit, "blblogreg")
  expect_type(sigma.blblogreg(fit), "double")
  expect_equal(length(coef.blblogreg(fit)), 4)
  expect_equal(dim(confint.blblogreg(fit)), c(3,2))
  expect_equal(length(predict.blblogreg(fit, new_data = data.frame("Lag1" = c(0.382, 0.4, 0.5),"Lag2" = c(1.03, 1.4, 0.5), "Volume" = c(1.2, 1.34, 1.1)))), 3)
})
