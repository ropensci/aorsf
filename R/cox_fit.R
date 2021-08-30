
cox_fit_center <- function(x,
                           y,
                           wts,
                           mtry,
                           n_vars_lc){

  fit <-
    survival::coxph.fit(
      x = x,
      y = y,
      strata = NULL,
      weights = wts,
      control = survival::coxph.control(toler.inf = Inf),
      method = 'efron',
      rownames = NULL,
      resid = FALSE,
      nocenter = c(-1,0,1)
    )

  coef_keep <- !is.na(fit$coefficients) & fit$coefficients != 0

  beta <- fit$coefficients[coef_keep]
  se <- sqrt(diag(fit$var)[coef_keep])
  pv <- pchisq((beta/se)^2, df = 1, lower.tail = FALSE)

  keep <- seq(min(length(pv),n_vars_lc))
  keep_cols <- match(keep, order(pv))

  fit <- survival::coxph.fit(
    x = x[, keep_cols, drop = FALSE],
    y = y,
    strata = NULL,
    weights = wts,
    control = survival::coxph.control(toler.inf = Inf),
    init = coef(fit)[keep_cols],
    method = 'efron',
    rownames = NULL,
    resid = FALSE,
    nocenter = c(-1,0,1)
  )

  list(
    coef = fit$coef,
    col_index = keep_cols
  )

}

cox_fit_nocenter <- function(x,
                             y,
                             wts,
                             mtry,
                             n_vars_lc){

  fit <-
    survival::coxph.fit(
      x = x,
      y = y,
      strata = NULL,
      weights = wts,
      control = survival::coxph.control(toler.inf = Inf),
      method = 'efron',
      rownames = NULL,
      resid = FALSE
    )

  if(mtry == n_vars_lc) return(
    list(
      coef = fit$coef,
      col_index = seq(ncol(x))
    )
  )

  coef_keep <- !is.na(fit$coefficients) & fit$coefficients != 0

  beta <- fit$coefficients[coef_keep]
  se <- sqrt(diag(fit$var)[coef_keep])
  pv <- pchisq((beta/se)^2, df = 1, lower.tail = FALSE)

  keep <- seq(min(length(pv),n_vars_lc))
  keep_cols <- match(keep, order(pv))

  fit <- survival::coxph.fit(
    x = x[, keep_cols, drop = FALSE],
    y = y,
    strata = NULL,
    weights = wts,
    control = survival::coxph.control(toler.inf = Inf),
    init = coef(fit)[keep_cols],
    method = 'efron',
    rownames = NULL,
    resid = FALSE
  )

  list(
    coef = fit$coef,
    lc = fit$linear.predictors+sum(fit$coef*fit$means),
    col_index = keep_cols
  )

}
