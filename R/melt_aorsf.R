
# need to make this to avoid possible memory leak in data.table melt

melt_aorsf <-
 function(data,
          id.vars,
          measure.vars,
          variable.name = "variable",
          value.name = "value") {

  # Select id variables and measure variables
  id_data <- select_cols(data, id.vars)

  measure_data <- select_cols(data, measure.vars)

  # Create a sequence variable to represent the variable names
  variable_data <- rep(names(measure_data), each = nrow(data))

  # Reshape the data
  long_data <- data.frame(id_data,
                          variable = variable_data,
                          value = unlist(measure_data, use.names = FALSE))

  names(long_data)[names(long_data) == 'variable'] <- variable.name
  names(long_data)[names(long_data) == 'value'] <- value.name

  return(long_data)
 }

