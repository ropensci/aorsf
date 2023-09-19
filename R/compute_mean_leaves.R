

compute_mean_leaves <- function(forest){

 if(is.null(forest$leaf_summary)){
  return(0)
 }

 collapse::fmean(
  vapply(forest$leaf_summary,
         function(leaf_smry) sum(leaf_smry != 0),
         FUN.VALUE = integer(1))
 )

}
