#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List call_coxph(NumericMatrix X, NumericVector time, NumericVector status) {
  // Load the survival package
  Function require("require");
  require("survival");

  // Create data frame from input
  int n = X.nrow(), p = X.ncol();
  List df(p + 2);
  CharacterVector col_names(p + 2);
  
  df[0] = time;
  col_names[0] = "time";
  df[1] = status;
  col_names[1] = "status";
  
  for (int j = 0; j < p; ++j) {
    df[2 + j] = X(_, j);
    col_names[2 + j] = "X" + std::to_string(j + 1);
  }
  
  df.attr("names") = col_names;
  df.attr("class") = "data.frame";
  df.attr("row.names") = seq(1, n);

  // Create formula: Surv(time, status) ~ .
  Function as_formula("as.formula");
  SEXP formula = as_formula("Surv(time, status) ~ .");

  // Call coxph
  Function coxph("coxph");
  List fit = coxph(_["formula"] = formula, _["data"] = df);
  
  return fit;
}