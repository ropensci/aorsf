
data <- survival::pbc

numeric_cols <- names(which(vapply(data, is.numeric, logical(1))))
categorical_cols <- names(which(vapply(data, is.factor, logical(1))))

impute_values <- c(get_means(data, numeric_cols),
                   get_modes(data, categorical_cols))

for(i in names(impute_values)){

 if(i %in% numeric_cols)
  expect_equal(impute_values[[i]], mean(data[[i]], na.rm = TRUE))

 if(i %in% categorical_cols)
  expect_equal(impute_values[[i]], mode(data[[i]]))

}

missing_index <- is.na(data)
colnames(missing_index) <- colnames(data)

data_imputed <- impute_meanmode(data,
                                cols = names(data),
                                values = impute_values)

test_that(
 desc = "means and modes are placed into data after imputation",
 code = {
  for(col in names(data)){

   if( any(missing_index[, col]) ){

    rows_imputed <- which(missing_index[, col])

    expect_equal(
     data_imputed[rows_imputed, col],
     rep(impute_values[[col]], length(rows_imputed))
    )
   }
  }
 }
)

