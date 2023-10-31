// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// coxph_fit_exported
List coxph_fit_exported(arma::mat& x_node, arma::mat& y_node, arma::vec& w_node, int method, double epsilon, arma::uword iter_max);
RcppExport SEXP _aorsf_coxph_fit_exported(SEXP x_nodeSEXP, SEXP y_nodeSEXP, SEXP w_nodeSEXP, SEXP methodSEXP, SEXP epsilonSEXP, SEXP iter_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x_node(x_nodeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y_node(y_nodeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w_node(w_nodeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type iter_max(iter_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(coxph_fit_exported(x_node, y_node, w_node, method, epsilon, iter_max));
    return rcpp_result_gen;
END_RCPP
}
// linreg_fit_exported
arma::mat linreg_fit_exported(arma::mat& x_node, arma::mat& y_node, arma::vec& w_node, bool do_scale, double epsilon, arma::uword iter_max);
RcppExport SEXP _aorsf_linreg_fit_exported(SEXP x_nodeSEXP, SEXP y_nodeSEXP, SEXP w_nodeSEXP, SEXP do_scaleSEXP, SEXP epsilonSEXP, SEXP iter_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x_node(x_nodeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y_node(y_nodeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w_node(w_nodeSEXP);
    Rcpp::traits::input_parameter< bool >::type do_scale(do_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type iter_max(iter_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(linreg_fit_exported(x_node, y_node, w_node, do_scale, epsilon, iter_max));
    return rcpp_result_gen;
END_RCPP
}
// logreg_fit_exported
arma::mat logreg_fit_exported(arma::mat& x_node, arma::mat& y_node, arma::vec& w_node, bool do_scale, double epsilon, arma::uword iter_max);
RcppExport SEXP _aorsf_logreg_fit_exported(SEXP x_nodeSEXP, SEXP y_nodeSEXP, SEXP w_nodeSEXP, SEXP do_scaleSEXP, SEXP epsilonSEXP, SEXP iter_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x_node(x_nodeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y_node(y_nodeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w_node(w_nodeSEXP);
    Rcpp::traits::input_parameter< bool >::type do_scale(do_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type iter_max(iter_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(logreg_fit_exported(x_node, y_node, w_node, do_scale, epsilon, iter_max));
    return rcpp_result_gen;
END_RCPP
}
// compute_cstat_surv_exported_vec
double compute_cstat_surv_exported_vec(arma::mat& y, arma::vec& w, arma::vec& p, bool pred_is_risklike);
RcppExport SEXP _aorsf_compute_cstat_surv_exported_vec(SEXP ySEXP, SEXP wSEXP, SEXP pSEXP, SEXP pred_is_risklikeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type pred_is_risklike(pred_is_risklikeSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_cstat_surv_exported_vec(y, w, p, pred_is_risklike));
    return rcpp_result_gen;
END_RCPP
}
// compute_cstat_surv_exported_uvec
double compute_cstat_surv_exported_uvec(arma::mat& y, arma::vec& w, arma::uvec& g, bool pred_is_risklike);
RcppExport SEXP _aorsf_compute_cstat_surv_exported_uvec(SEXP ySEXP, SEXP wSEXP, SEXP gSEXP, SEXP pred_is_risklikeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type g(gSEXP);
    Rcpp::traits::input_parameter< bool >::type pred_is_risklike(pred_is_risklikeSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_cstat_surv_exported_uvec(y, w, g, pred_is_risklike));
    return rcpp_result_gen;
END_RCPP
}
// compute_cstat_clsf_exported
double compute_cstat_clsf_exported(arma::vec& y, arma::vec& w, arma::vec& p);
RcppExport SEXP _aorsf_compute_cstat_clsf_exported(SEXP ySEXP, SEXP wSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_cstat_clsf_exported(y, w, p));
    return rcpp_result_gen;
END_RCPP
}
// compute_logrank_exported
double compute_logrank_exported(arma::mat& y, arma::vec& w, arma::uvec& g);
RcppExport SEXP _aorsf_compute_logrank_exported(SEXP ySEXP, SEXP wSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_logrank_exported(y, w, g));
    return rcpp_result_gen;
END_RCPP
}
// compute_gini_exported
double compute_gini_exported(arma::mat& y, arma::vec& w, arma::uvec& g);
RcppExport SEXP _aorsf_compute_gini_exported(SEXP ySEXP, SEXP wSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_gini_exported(y, w, g));
    return rcpp_result_gen;
END_RCPP
}
// compute_pred_prob_exported
arma::vec compute_pred_prob_exported(arma::mat& y, arma::vec& w);
RcppExport SEXP _aorsf_compute_pred_prob_exported(SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_pred_prob_exported(y, w));
    return rcpp_result_gen;
END_RCPP
}
// expand_y_clsf
arma::mat expand_y_clsf(arma::vec& y, arma::uword n_class);
RcppExport SEXP _aorsf_expand_y_clsf(SEXP ySEXP, SEXP n_classSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_class(n_classSEXP);
    rcpp_result_gen = Rcpp::wrap(expand_y_clsf(y, n_class));
    return rcpp_result_gen;
END_RCPP
}
// is_col_splittable_exported
bool is_col_splittable_exported(arma::mat& x, arma::mat& y, arma::uvec& r, arma::uword j);
RcppExport SEXP _aorsf_is_col_splittable_exported(SEXP xSEXP, SEXP ySEXP, SEXP rSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(is_col_splittable_exported(x, y, r, j));
    return rcpp_result_gen;
END_RCPP
}
// find_cuts_survival_exported
List find_cuts_survival_exported(arma::mat& y, arma::vec& w, arma::vec& lincomb, double leaf_min_events, double leaf_min_obs, int split_rule_R);
RcppExport SEXP _aorsf_find_cuts_survival_exported(SEXP ySEXP, SEXP wSEXP, SEXP lincombSEXP, SEXP leaf_min_eventsSEXP, SEXP leaf_min_obsSEXP, SEXP split_rule_RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lincomb(lincombSEXP);
    Rcpp::traits::input_parameter< double >::type leaf_min_events(leaf_min_eventsSEXP);
    Rcpp::traits::input_parameter< double >::type leaf_min_obs(leaf_min_obsSEXP);
    Rcpp::traits::input_parameter< int >::type split_rule_R(split_rule_RSEXP);
    rcpp_result_gen = Rcpp::wrap(find_cuts_survival_exported(y, w, lincomb, leaf_min_events, leaf_min_obs, split_rule_R));
    return rcpp_result_gen;
END_RCPP
}
// sprout_node_survival_exported
List sprout_node_survival_exported(arma::mat& y, arma::vec& w);
RcppExport SEXP _aorsf_sprout_node_survival_exported(SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(sprout_node_survival_exported(y, w));
    return rcpp_result_gen;
END_RCPP
}
// find_rows_inbag_exported
arma::uvec find_rows_inbag_exported(arma::uvec rows_oobag, arma::uword n_obs);
RcppExport SEXP _aorsf_find_rows_inbag_exported(SEXP rows_oobagSEXP, SEXP n_obsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type rows_oobag(rows_oobagSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_obs(n_obsSEXP);
    rcpp_result_gen = Rcpp::wrap(find_rows_inbag_exported(rows_oobag, n_obs));
    return rcpp_result_gen;
END_RCPP
}
// x_submat_mult_beta_exported
arma::vec x_submat_mult_beta_exported(arma::mat& x, arma::mat& y, arma::vec& w, arma::uvec& x_rows, arma::uvec& x_cols, arma::vec& beta);
RcppExport SEXP _aorsf_x_submat_mult_beta_exported(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP x_rowsSEXP, SEXP x_colsSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type x_rows(x_rowsSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type x_cols(x_colsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(x_submat_mult_beta_exported(x, y, w, x_rows, x_cols, beta));
    return rcpp_result_gen;
END_RCPP
}
// x_submat_mult_beta_pd_exported
arma::vec x_submat_mult_beta_pd_exported(arma::mat& x, arma::mat& y, arma::vec& w, arma::uvec& x_rows, arma::uvec& x_cols, arma::vec& beta, arma::vec& pd_x_vals, arma::uvec& pd_x_cols);
RcppExport SEXP _aorsf_x_submat_mult_beta_pd_exported(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP x_rowsSEXP, SEXP x_colsSEXP, SEXP betaSEXP, SEXP pd_x_valsSEXP, SEXP pd_x_colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type x_rows(x_rowsSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type x_cols(x_colsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type pd_x_vals(pd_x_valsSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type pd_x_cols(pd_x_colsSEXP);
    rcpp_result_gen = Rcpp::wrap(x_submat_mult_beta_pd_exported(x, y, w, x_rows, x_cols, beta, pd_x_vals, pd_x_cols));
    return rcpp_result_gen;
END_RCPP
}
// scale_x_exported
List scale_x_exported(arma::mat& x, arma::vec& w);
RcppExport SEXP _aorsf_scale_x_exported(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(scale_x_exported(x, w));
    return rcpp_result_gen;
END_RCPP
}
// cph_scale
List cph_scale(arma::mat& x, arma::vec& w);
RcppExport SEXP _aorsf_cph_scale(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(cph_scale(x, w));
    return rcpp_result_gen;
END_RCPP
}
// orsf_cpp
List orsf_cpp(arma::mat& x, arma::mat& y, arma::vec& w, arma::uword tree_type_R, Rcpp::IntegerVector& tree_seeds, Rcpp::List& loaded_forest, Rcpp::RObject lincomb_R_function, Rcpp::RObject oobag_R_function, arma::uword n_tree, arma::uword mtry, bool sample_with_replacement, double sample_fraction, arma::uword vi_type_R, double vi_max_pvalue, double leaf_min_events, double leaf_min_obs, arma::uword split_rule_R, double split_min_events, double split_min_obs, double split_min_stat, arma::uword split_max_cuts, arma::uword split_max_retry, arma::uword lincomb_type_R, double lincomb_eps, arma::uword lincomb_iter_max, bool lincomb_scale, double lincomb_alpha, arma::uword lincomb_df_target, arma::uword lincomb_ties_method, bool pred_mode, arma::uword pred_type_R, arma::vec pred_horizon, bool pred_aggregate, bool oobag, arma::uword oobag_eval_type_R, arma::uword oobag_eval_every, int pd_type_R, std::vector<arma::mat>& pd_x_vals, std::vector<arma::uvec>& pd_x_cols, arma::vec& pd_probs, unsigned int n_thread, bool write_forest, bool run_forest, int verbosity);
RcppExport SEXP _aorsf_orsf_cpp(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP tree_type_RSEXP, SEXP tree_seedsSEXP, SEXP loaded_forestSEXP, SEXP lincomb_R_functionSEXP, SEXP oobag_R_functionSEXP, SEXP n_treeSEXP, SEXP mtrySEXP, SEXP sample_with_replacementSEXP, SEXP sample_fractionSEXP, SEXP vi_type_RSEXP, SEXP vi_max_pvalueSEXP, SEXP leaf_min_eventsSEXP, SEXP leaf_min_obsSEXP, SEXP split_rule_RSEXP, SEXP split_min_eventsSEXP, SEXP split_min_obsSEXP, SEXP split_min_statSEXP, SEXP split_max_cutsSEXP, SEXP split_max_retrySEXP, SEXP lincomb_type_RSEXP, SEXP lincomb_epsSEXP, SEXP lincomb_iter_maxSEXP, SEXP lincomb_scaleSEXP, SEXP lincomb_alphaSEXP, SEXP lincomb_df_targetSEXP, SEXP lincomb_ties_methodSEXP, SEXP pred_modeSEXP, SEXP pred_type_RSEXP, SEXP pred_horizonSEXP, SEXP pred_aggregateSEXP, SEXP oobagSEXP, SEXP oobag_eval_type_RSEXP, SEXP oobag_eval_everySEXP, SEXP pd_type_RSEXP, SEXP pd_x_valsSEXP, SEXP pd_x_colsSEXP, SEXP pd_probsSEXP, SEXP n_threadSEXP, SEXP write_forestSEXP, SEXP run_forestSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type tree_type_R(tree_type_RSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type tree_seeds(tree_seedsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type loaded_forest(loaded_forestSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type lincomb_R_function(lincomb_R_functionSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type oobag_R_function(oobag_R_functionSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_tree(n_treeSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type mtry(mtrySEXP);
    Rcpp::traits::input_parameter< bool >::type sample_with_replacement(sample_with_replacementSEXP);
    Rcpp::traits::input_parameter< double >::type sample_fraction(sample_fractionSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type vi_type_R(vi_type_RSEXP);
    Rcpp::traits::input_parameter< double >::type vi_max_pvalue(vi_max_pvalueSEXP);
    Rcpp::traits::input_parameter< double >::type leaf_min_events(leaf_min_eventsSEXP);
    Rcpp::traits::input_parameter< double >::type leaf_min_obs(leaf_min_obsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type split_rule_R(split_rule_RSEXP);
    Rcpp::traits::input_parameter< double >::type split_min_events(split_min_eventsSEXP);
    Rcpp::traits::input_parameter< double >::type split_min_obs(split_min_obsSEXP);
    Rcpp::traits::input_parameter< double >::type split_min_stat(split_min_statSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type split_max_cuts(split_max_cutsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type split_max_retry(split_max_retrySEXP);
    Rcpp::traits::input_parameter< arma::uword >::type lincomb_type_R(lincomb_type_RSEXP);
    Rcpp::traits::input_parameter< double >::type lincomb_eps(lincomb_epsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type lincomb_iter_max(lincomb_iter_maxSEXP);
    Rcpp::traits::input_parameter< bool >::type lincomb_scale(lincomb_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type lincomb_alpha(lincomb_alphaSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type lincomb_df_target(lincomb_df_targetSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type lincomb_ties_method(lincomb_ties_methodSEXP);
    Rcpp::traits::input_parameter< bool >::type pred_mode(pred_modeSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type pred_type_R(pred_type_RSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pred_horizon(pred_horizonSEXP);
    Rcpp::traits::input_parameter< bool >::type pred_aggregate(pred_aggregateSEXP);
    Rcpp::traits::input_parameter< bool >::type oobag(oobagSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type oobag_eval_type_R(oobag_eval_type_RSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type oobag_eval_every(oobag_eval_everySEXP);
    Rcpp::traits::input_parameter< int >::type pd_type_R(pd_type_RSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat>& >::type pd_x_vals(pd_x_valsSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::uvec>& >::type pd_x_cols(pd_x_colsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type pd_probs(pd_probsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_thread(n_threadSEXP);
    Rcpp::traits::input_parameter< bool >::type write_forest(write_forestSEXP);
    Rcpp::traits::input_parameter< bool >::type run_forest(run_forestSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(orsf_cpp(x, y, w, tree_type_R, tree_seeds, loaded_forest, lincomb_R_function, oobag_R_function, n_tree, mtry, sample_with_replacement, sample_fraction, vi_type_R, vi_max_pvalue, leaf_min_events, leaf_min_obs, split_rule_R, split_min_events, split_min_obs, split_min_stat, split_max_cuts, split_max_retry, lincomb_type_R, lincomb_eps, lincomb_iter_max, lincomb_scale, lincomb_alpha, lincomb_df_target, lincomb_ties_method, pred_mode, pred_type_R, pred_horizon, pred_aggregate, oobag, oobag_eval_type_R, oobag_eval_every, pd_type_R, pd_x_vals, pd_x_cols, pd_probs, n_thread, write_forest, run_forest, verbosity));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_aorsf_coxph_fit_exported", (DL_FUNC) &_aorsf_coxph_fit_exported, 6},
    {"_aorsf_linreg_fit_exported", (DL_FUNC) &_aorsf_linreg_fit_exported, 6},
    {"_aorsf_logreg_fit_exported", (DL_FUNC) &_aorsf_logreg_fit_exported, 6},
    {"_aorsf_compute_cstat_surv_exported_vec", (DL_FUNC) &_aorsf_compute_cstat_surv_exported_vec, 4},
    {"_aorsf_compute_cstat_surv_exported_uvec", (DL_FUNC) &_aorsf_compute_cstat_surv_exported_uvec, 4},
    {"_aorsf_compute_cstat_clsf_exported", (DL_FUNC) &_aorsf_compute_cstat_clsf_exported, 3},
    {"_aorsf_compute_logrank_exported", (DL_FUNC) &_aorsf_compute_logrank_exported, 3},
    {"_aorsf_compute_gini_exported", (DL_FUNC) &_aorsf_compute_gini_exported, 3},
    {"_aorsf_compute_pred_prob_exported", (DL_FUNC) &_aorsf_compute_pred_prob_exported, 2},
    {"_aorsf_expand_y_clsf", (DL_FUNC) &_aorsf_expand_y_clsf, 2},
    {"_aorsf_is_col_splittable_exported", (DL_FUNC) &_aorsf_is_col_splittable_exported, 4},
    {"_aorsf_find_cuts_survival_exported", (DL_FUNC) &_aorsf_find_cuts_survival_exported, 6},
    {"_aorsf_sprout_node_survival_exported", (DL_FUNC) &_aorsf_sprout_node_survival_exported, 2},
    {"_aorsf_find_rows_inbag_exported", (DL_FUNC) &_aorsf_find_rows_inbag_exported, 2},
    {"_aorsf_x_submat_mult_beta_exported", (DL_FUNC) &_aorsf_x_submat_mult_beta_exported, 6},
    {"_aorsf_x_submat_mult_beta_pd_exported", (DL_FUNC) &_aorsf_x_submat_mult_beta_pd_exported, 8},
    {"_aorsf_scale_x_exported", (DL_FUNC) &_aorsf_scale_x_exported, 2},
    {"_aorsf_cph_scale", (DL_FUNC) &_aorsf_cph_scale, 2},
    {"_aorsf_orsf_cpp", (DL_FUNC) &_aorsf_orsf_cpp, 44},
    {NULL, NULL, 0}
};

RcppExport void R_init_aorsf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
