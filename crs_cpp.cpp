#include <RcppArmadillo.h>
#include <utility>
#include <tuple>
#include <string>

using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//' Armadillo-based reimplementation of `stats::quantile()` for `type=7`.
//' @description Allows for NAs.
//' @param x a vec of values for which to find the quantiles
//' @param probs a vec of probilities that specify which quantiles 
//' @returns a vec of quantiles
//' @examples
//' set.seed(12345)
//' y=rnorm(1e5)
//' quantile(y)-quantile_type7(y,seq(0,1,length.out=5))
//' quantile(y, 0.975) - quantile_type7(y, 0.975)
//' quantile(y,  probs = c(0.1, 0.5, 1, 2, 5, 10, 50, NA)/100) -
//'     quantile_type7(y,  probs = c(0.1, 0.5, 1, 2, 5, 10, 50, NA)/100)
// [[Rcpp::export]]
vec quantile_type7(vec x, vec probs) { 
  x = sort(x);
  size_t n = x.size();
  size_t m = probs.size();
  vec q(m);
  if (n == 0) {
    for (size_t i=0; i<m; i++)
      q[i] = arma::datum::nan;
    return q;
  } else 
    for (size_t i=0; i<m; i++) {
      if (std::isnan(probs[i])) {
	q[i] = arma::datum::nan;
      } else {
	double index = 1 + (n - 1) * probs[i];
	int lo = std::floor(index) - 1; // 0-based
	int hi = std::ceil(index) - 1;  // 0-based
	double gamma = index - std::floor(index);
	double q_lo = x[lo];
	double q_hi = (hi > lo ? x[hi] : q_lo);
	q[i] = (1 - gamma) * q_lo + gamma * q_hi;
      }
    }
  return q;
}

mat chol2inv(mat L) {
  mat L_inv = inv(L);
  mat A_inv = L_inv.t() * L_inv;
  return A_inv;
}

// penalty matrix
std::pair<mat,mat> crs_FP_cpp(vec knots, vec h) {
  size_t k = knots.size();
  mat B(k-2, k-2), D(k-2,k);
  for (size_t i=0; i<(k-2); i++) {
    D(i,i) = 1.0/h(i);
    D(i,i+1) = -1.0/h(i)-1.0/h(i+1);
    D(i,i+2) = 1.0/h(i+1);
    B(i,i) = (h(i)+h(i+1))/3.0;
    if (i < (k-2-1)) { B(i,i+1) = h(i+1)/6.0; B(i+1,i) = h(i+1)/6.0; }
  }
  mat F_mat = chol2inv(chol(B,"lower")) * D;
  mat P_mat = D.t() * F_mat;
  P_mat = (P_mat+P_mat.t())*0.5;
  return std::make_pair(F_mat,P_mat);
}

// // Bug when derivs==1
// std::pair<mat,mat> crs_cpp2(vec x, vec knots, size_t derivs=0, bool intercept=true) {
//   size_t n = x.size();
//   unsigned int k = knots.size();
//   vec h = diff(knots);
//   std::pair<mat,mat> F_P = crs_FP_cpp(knots,h);
//   mat F_mat = F_P.first;
//   mat F_mat1 = join_cols(zeros(k).t(),F_mat);
//   mat F_mat2 = join_cols(F_mat,zeros(k).t());
//   mat P_mat = F_P.second;
//   uvec condition_min, condition_max;
//   vec x_min, x_max;
//   if (derivs==0) {
//     condition_min = find(x < min(knots));
//     condition_max = find(x > max(knots));
//     x_min = x(condition_min);
//     x_max = x(condition_max);
//     x(condition_min) = x(condition_min)*0.0 + min(knots);
//     x(condition_max) = x(condition_max)*0.0 + max(knots);
//   } else {
//     condition_min = find(x <= min(knots));
//     condition_max = find(x >= max(knots));
//   }
//   mat condition(n, k-1), a_minus(n, k-1), a_plus(n, k-1), c_minus(n, k-1), c_plus(n, k-1);
//   for(unsigned int l=0; l<(k-1); l++) {
//     uvec col{l};
//     condition(find(x>=knots(l) && x<knots(l+1)), col) += 1.0;
//     if (derivs==0) {
//       a_minus.col(l)=(knots(l+1)-x)/h(l);
//       a_plus.col(l) = (x-knots(l))/h(l);
//       c_minus.col(l) = (pow(knots(l+1)-x,3)/h(l) - h(l)*(knots(l+1)-x))/6.0;
//       c_plus.col(l) = (pow(x-knots(l),3)/h(l) - h(l)*(x-knots(l)))/6.0;
//     } else if (derivs==1) {
//       a_minus.col(l) -= 1/h(l);
//       a_plus.col(l)  += 1/h(l);
//       c_minus.col(l) = (-3*pow(knots(l+1)-x,2)/h(l) + h(l))/6.0;
//       c_plus.col(l)  = (3*pow(x-knots(l),2)/h(l) - h(l))/6.0;
//     } else {
//       c_minus.col(l) = (knots(l+1)-x)/h(l);
//       c_plus.col(l)  = (x-knots(l))/h(l);
//     }
//   }
//   a_minus %= condition;
//   a_plus %= condition;
//   c_minus %= condition;
//   c_plus %= condition;
//   mat Ident = eye(k-1,k-1);
//   mat Mat_j = join_rows(Ident, zeros(k-1));
//   mat Mat_j_1 = join_rows(zeros(k-1), Ident);
//   mat b = c_minus * F_mat1 + c_plus*F_mat2 + a_minus * Mat_j + a_plus*Mat_j_1;
//   if (derivs==0) {
//     uvec x_max_knots = find(x == max(knots));
//     if (x_max_knots.size()>0) {
//       uvec col{k-1};
//       b(x_max_knots,col) = b(x_max_knots,col)*0.0+1.0;
//     }
//   }
//   if (condition_min.size()>0) {
//     vec v1 = x_min - min(knots);
//     vec e0=zeros(k); e0(0)=1.0;
//     vec e1=zeros(k); e1(1)=1.0;
//     rowvec v2 = -h(0)/6.0 * F_mat.row(0)-1.0/h(0)*e0.t() + 1.0/h(0)*e1.t();
//     for (size_t i=0; i<condition_min.size(); i++)
//       b.row(condition_min(i)) += (derivs==0) ? v1(i)*v2 :
// 	(derivs==1) ? v2 : zeros(k).t();
//   }
//   if (condition_max.size()>0) {
//     vec v1 = x_max - max(knots);
//     vec ek1=zeros(k); ek1(k-1)=1.0;
//     vec ek2=zeros(k); ek2(k-2)=1.0;
//     rowvec v2 = h(k-2)/6.0 * F_mat.row(k-3)-1.0/h(k-2)*ek2.t() + 1.0/h(k-2)*ek1.t();
//     for (size_t i=0; i<condition_max.size(); i++)
//       b.row(condition_max(i)) += (derivs==0) ? v1(i)*v2 :
// 	(derivs==1) ? v2 : zeros(k).t();
//   }
//   if (!intercept) {
//     b = b(span::all, span(1,b.n_cols-1));
//     P_mat = P_mat(span(1, P_mat.n_rows-1), span(1, P_mat.n_cols-1));
//   }
//   return std::make_pair(b, P_mat);
// }

// function evaluation called in crs_cpp
std::pair<mat,mat> crs_cpp_D0(vec x, vec knots, bool intercept=true) {
  size_t n = x.size();
  unsigned int k = knots.size();
  vec h = diff(knots);
  std::pair<mat,mat> F_P = crs_FP_cpp(knots,h);
  mat F_mat = F_P.first;
  mat F_mat1 = join_cols(zeros(k).t(),F_mat);
  mat F_mat2 = join_cols(F_mat,zeros(k).t());
  mat P_mat = F_P.second;
  uvec condition_min = find(x < min(knots));
  uvec condition_max = find(x > max(knots));
  vec x_min = x(condition_min);
  vec x_max = x(condition_max);
  x(condition_min) = x(condition_min)*0.0 + min(knots);
  x(condition_max) = x(condition_max)*0.0 + max(knots);
  mat condition(n, k-1), a_minus(n, k-1), a_plus(n, k-1), c_minus(n, k-1), c_plus(n, k-1);
  for(unsigned int l=0; l<(k-1); l++) {
    uvec col{l};
    condition(find(x>=knots(l) && x<knots(l+1)), col) += 1.0;
    a_minus.col(l)=(knots(l+1)-x)/h(l);
    a_plus.col(l) = (x-knots(l))/h(l);
    c_minus.col(l) = (pow(knots(l+1)-x,3)/h(l) - h(l)*(knots(l+1)-x))/6.0;
    c_plus.col(l) = (pow(x-knots(l),3)/h(l) - h(l)*(x-knots(l)))/6.0;
  }
  a_minus %= condition;
  a_plus %= condition;
  c_minus %= condition;
  c_plus %= condition;
  mat Ident = eye(k-1,k-1);
  mat Mat_j = join_rows(Ident, zeros(k-1));
  mat Mat_j_1 = join_rows(zeros(k-1), Ident);
  mat b = c_minus * F_mat1 + c_plus*F_mat2 + a_minus * Mat_j + a_plus*Mat_j_1;
  uvec x_max_knots = find(x == max(knots));
  if (x_max_knots.size()>0) {
    uvec col{k-1};
    b(x_max_knots,col) = b(x_max_knots,col)*0.0+1.0;
  }
  if (condition_min.size()>0) {
    vec v1 = x_min - min(knots);
    vec e0=zeros(k); e0(0)=1.0;
    vec e1=zeros(k); e1(1)=1.0;
    rowvec v2 = -h(0)/6.0 * F_mat.row(0)-1.0/h(0)*e0.t() + 1.0/h(0)*e1.t();
    for (size_t i=0; i<condition_min.size(); i++) b.row(condition_min(i)) += v1(i)*v2;
  }
  if (condition_max.size()>0) {
    vec v1 = x_max - max(knots);
    vec ek1=zeros(k); ek1(k-1)=1.0;
    vec ek2=zeros(k); ek2(k-2)=1.0;
    rowvec v2 = h(k-2)/6.0 * F_mat.row(k-3)-1.0/h(k-2)*ek2.t() + 1.0/h(k-2)*ek1.t();
    for (size_t i=0; i<condition_max.size(); i++) b.row(condition_max(i)) += v1(i)*v2;
  }
  if (!intercept) {
    b = b(span::all, span(1,b.n_cols-1));
    P_mat = P_mat(span(1, P_mat.n_rows-1), span(1, P_mat.n_cols-1));
  }
  return std::make_pair(b, P_mat);
}

// first derivative called in crs_cpp
std::pair<mat,mat> crs_cpp_D1(vec x, vec knots, bool intercept=true) {
  size_t n = x.size();
  unsigned int k = knots.size();
  vec h = diff(knots);
  std::pair<mat,mat> F_P = crs_FP_cpp(knots,h);
  mat F_mat = F_P.first;
  mat F_mat1 = join_cols(zeros(k).t(),F_mat);
  mat F_mat2 = join_cols(F_mat,zeros(k).t());
  mat P_mat = F_P.second;
  uvec condition_min = find(x < min(knots));
  uvec condition_max = find(x > max(knots));
  mat condition(n, k-1), a_minus(n, k-1), a_plus(n, k-1), c_minus(n, k-1), c_plus(n, k-1);
  for(unsigned int l=0; l<(k-1); l++) {
    uvec col{l};
    condition(find(x>=knots(l) && x<knots(l+1)), col) += 1.0;
    // a_minus.col(l)=(knots(l+1)-x)/h(l);
    // a_plus.col(l) = (x-knots(l))/h(l);
    // c_minus.col(l) = (pow(knots(l+1)-x,3)/h(l) - h(l)*(knots(l+1)-x))/6.0;
    // c_plus.col(l) = (pow(x-knots(l),3)/h(l) - h(l)*(x-knots(l)))/6.0;
    a_minus.col(l) -= 1/h(l);
    a_plus.col(l)  += 1/h(l);
    c_minus.col(l) = (-3*pow(knots(l+1)-x,2)/h(l) + h(l))/6.0;
    c_plus.col(l)  = (3*pow(x-knots(l),2)/h(l) - h(l))/6.0;
  }
  a_minus %= condition;
  a_plus %= condition;
  c_minus %= condition;
  c_plus %= condition;
  mat Ident = eye(k-1,k-1);
  mat Mat_j = join_rows(Ident, zeros(k-1));
  mat Mat_j_1 = join_rows(zeros(k-1), Ident);
  mat b = c_minus * F_mat1 + c_plus*F_mat2 + a_minus * Mat_j + a_plus*Mat_j_1;
  if (condition_min.size()>0) {
    vec e0=zeros(k); e0(0)=1.0;
    vec e1=zeros(k); e1(1)=1.0;
    rowvec v2 = -h(0)/6.0 * F_mat.row(0)-1.0/h(0)*e0.t() + 1.0/h(0)*e1.t();
    for (size_t i=0; i<condition_min.size(); i++) b.row(condition_min(i)) = v2;
  }
  if (condition_max.size()>0) {
    vec ek1=zeros(k); ek1(k-1)=1.0;
    vec ek2=zeros(k); ek2(k-2)=1.0;
    rowvec v2 = h(k-2)/6.0 * F_mat.row(k-3)-1.0/h(k-2)*ek2.t() + 1.0/h(k-2)*ek1.t();
    for (size_t i=0; i<condition_max.size(); i++) b.row(condition_max(i)) = v2;
  }
  if (!intercept) {
    b = b(span::all, span(1,b.n_cols-1));
    P_mat = P_mat(span(1, P_mat.n_rows-1), span(1, P_mat.n_cols-1));
  }
  return std::make_pair(b, P_mat);
}

// second derivative called in crs_cpp
std::pair<mat,mat> crs_cpp_D2(vec x, vec knots, bool intercept=true) {
  size_t n = x.size();
  unsigned int k = knots.size();
  vec h = diff(knots);
  std::pair<mat,mat> F_P = crs_FP_cpp(knots,h);
  mat F_mat = F_P.first;
  mat F_mat1 = join_cols(zeros(k).t(),F_mat);
  mat F_mat2 = join_cols(F_mat,zeros(k).t());
  mat P_mat = F_P.second;
  uvec condition_min = find(x < min(knots));
  uvec condition_max = find(x > max(knots));
  mat condition(n, k-1), a_minus(n, k-1), a_plus(n, k-1), c_minus(n, k-1), c_plus(n, k-1);
  for(unsigned int l=0; l<(k-1); l++) {
    uvec col{l};
    condition(find(x>=knots(l) && x<knots(l+1)), col) += 1.0;
    // a_minus.col(l)=(knots(l+1)-x)/h(l);
    // a_plus.col(l) = (x-knots(l))/h(l);
    // c_minus.col(l) = (pow(knots(l+1)-x,3)/h(l) - h(l)*(knots(l+1)-x))/6.0;
    // c_plus.col(l) = (pow(x-knots(l),3)/h(l) - h(l)*(x-knots(l)))/6.0;
    c_minus.col(l) = (knots(l+1)-x)/h(l);
    c_plus.col(l)  = (x-knots(l))/h(l);
  }
  a_minus %= condition;
  a_plus %= condition;
  c_minus %= condition;
  c_plus %= condition;
  mat Ident = eye(k-1,k-1);
  mat Mat_j = join_rows(Ident, zeros(k-1));
  mat Mat_j_1 = join_rows(zeros(k-1), Ident);
  mat b = c_minus * F_mat1 + c_plus*F_mat2 + a_minus * Mat_j + a_plus*Mat_j_1;
  if (condition_min.size()>0) {
    vec e0=zeros(k);
    for (size_t i=0; i<condition_min.size(); i++) b.row(condition_min(i)) = e0.t();
  }
  if (condition_max.size()>0) {
    vec e0=zeros(k);
    for (size_t i=0; i<condition_max.size(); i++) b.row(condition_max(i)) = e0.t();
  }
  if (!intercept) {
    b = b(span::all, span(1,b.n_cols-1));
    P_mat = P_mat(span(1, P_mat.n_rows-1), span(1, P_mat.n_cols-1));
  }
  return std::make_pair(b, P_mat);
}

//' C++ reimplementation of survPen::crs, including derivatives.
//' @description The re-implementation uses the same algorithm but uses
//' different defaults, and allows for derivatives.
//' For the changes in defaults: the function here assumes that
//' the intercept is false (as per `splines::ns`); when the knots are
//' not specified, the function defaults to using `x` rather than `unique(x)`
//' for the quantiles; and the `df` argument comes before the `knots` argument.
//' To use the `survPen::crs` defaults in R, use
//' `crs_cpp(..., intercept=TRUE, survPen_compatible=TRUE)`.
//' @param x a vec of values to provide
//' @param df a non-negative integer for the degrees of freedom (ignored if
//' the knots are specified)
//' @param knots a vec of knots (defaults to NULL)
//' @param derivs a non-negative integer for the derivative (default 0)
//' @param intercept a bool for whether to include the intercept
//' @param survPen_compatible a bool for whether to use `unique(x)` for determining the knots if the `knots` are not specified 
//' @return an `Rcpp::NumericMatrix` with attributes `class="crs_cpp"`,
//' `pen` for the penalty matrix, `knots` for the knots, `derivs` for the order of
//' derivatives, and `intercept` for a bool/logical for whether the intercept is
//' included.
//' @examples
//' library(survPen)
//' x = seq(0,1,l=11)
//' knots = c(1,3,5,7)/10
//' crs1 = crs_cpp(x,knots=knots)
//' crs2 = crs(x,knots,intercept=FALSE)
//' range(crs1 - crs2$bs)
//' range(attr(crs1,"pen") - crs2$pen)
//' all(predict(crs1, rev(x)) == crs1[length(x):1,])
//' crs3 = crs_cpp(x,knots=knots,intercept=TRUE)
//' crs4 = crs(x,knots)
//' range(crs3 - crs4$bs)
//' range(attr(crs3,"pen") - crs4$pen)
//' all(update(crs1,x,intercept=TRUE) == crs3)
// [[Rcpp::export]]
Rcpp::NumericMatrix crs_cpp(vec x,
			    size_t df = 10,
			    Rcpp::Nullable<vec> knots = R_NilValue,
			    size_t derivs = 0,
			    bool intercept=false,
			    bool survPen_compatible = false) {
  size_t n = x.size();
  vec knotsi;
  if (knots.isNull()) {
    if (df<3)
      Rcpp::stop("Number of knots should be at least 3, 1 interior plus 2 boundaries");
    if (n<2)
      Rcpp::stop("Please specify at least 2 values or specify at least 3 knots via knots=...");
    vec probs = linspace(0.0,1.0,df);
    vec x_copy = survPen_compatible ? unique(x) : x;
    knotsi = quantile_type7(x_copy, probs);
  } else {
    knotsi = Rcpp::as<arma::vec>(knots);
    df = knotsi.size();
    if (df<3)
      Rcpp::stop("Please specify at least 3 knots, 1 interior plus 2 boundaries");
  }
  std::pair<mat,mat> pr =
    (derivs==0) ? crs_cpp_D0(x,knotsi,intercept) :
    (derivs==1) ? crs_cpp_D1(x,knotsi,intercept) :
    crs_cpp_D2(x,knotsi,intercept);
  Rcpp::NumericMatrix X = Rcpp::wrap(pr.first);
  X.attr("class") = Rcpp::CharacterVector{"crs_cpp","matrix"};
  X.attr("pen") = pr.second;
  X.attr("knots") = knotsi;
  X.attr("derivs") = derivs;
  X.attr("intercept") = intercept;
  return X;
}

//' C++ reimplementation of `smf`, including derivatives.
//' @description The re-implementation uses the same algorithm but uses
//' different defaults. Note that the inclusion of an
//' intercept is the same as `crs_cpp`, while no intercept uses a projection
//' based on a QR decomposition.
//' For the changes in defaults: the function here assumes that
//' the intercept is false (as per `splines::ns`); and when the knots are
//' not specified, the function defaults to using `x` rather than `unique(x)`
//' for the quantiles. The `survPen` implementation for `smf` assumes no
//' intercept unless there is a `by` variable. A `by` variable `x2` would
//' be implemented here by an interaction, such as
//' `smf_cpp(x,knots,intercept=TRUE):x2`; see the examples. 
//' @param x a vec of values to provide
//' @param df a non-negative integer for the degrees of freedom (ignored if
//' the knots are specified)
//' @param knots a vec of knots (defaults to NULL)
//' @param derivs a non-negative integer for the derivative (default 0)
//' @param intercept a bool for whether to include the intercept
//' @param survPen_compatible a bool for whether to use `unique(x)` for
//' determining the knots if the `knots` are not specified 
//' @return an `Rcpp::NumericMatrix` with attributes `class="smf_cpp"`,
//' `pen` for the penalty matrix, `knots` for the knots, `derivs` for the order of
//' derivatives, and `intercept` for a bool/logical for whether the intercept is
//' included. If there is no intercept, then there will also be a `Z` attribute
//' which is the QR-based projection matrix.
//' @suggests survPen
//' @examples
//' @rd smf_cpp
//' ## show that we get the same as survPen (within numerical precision)
//' library(survPen)
//' x = seq(0,1,l=11)
//' knots = c(1,3,5,7)/10
//' smf1 = smf_cpp(x,knots=knots)
//' smf2 = smooth.cons("x",knots=list(knots),df=4,option="smf",
//'                    data.spec=data.frame(x),name="smf(time)")
//' smf1 - smf2$X
//' range(smf1 - smf2$X)
//' range(attr(smf1,"pen") - smf2$pen[[1]])
//' range(attr(smf1, "Z") - smf2$Z.smf[[1]])
//' all(predict(smf1, rev(x)) == smf1[length(x):1,])
//' 
//' x2 = x*0+2
//' smf3 = smf_cpp(x,knots=knots,intercept=TRUE)*x2
//' smf4 = smooth.cons("x",knots=list(knots),df=4,option="smf",by="x2",
//'                    data.spec=data.frame(x,x2),name="smf(time)")
//' range(smf3 - smf4$X)
//' range(attr(smf3,"pen") - smf4$pen[[1]])
//' ## smf3 does not have a "Z" attribute
//' is.null(attr(smf3,"Z"))
//' all(smf3 == crs_cpp(x,knots=knots,intercept=TRUE)*x2)
//' all(model.matrix(~smf_cpp(x,knots=knots,intercept=TRUE):x2-1,
//'                  newdata=data.frame(x,x2)) == smf3)
//' all(update(smf1,x,intercept=TRUE)*x2 == smf3)
//' 
//' ## Check derivatives using finite differences
//' smfD = function(t,...,eps=1e-5) {
//'     object = smf_cpp(t,...)
//'     (predict(object, t+eps) - predict(object, t-eps))/(2*eps)
//' }
//' smfD2 = function(t,...,eps=1e-5) {
//'     object = smf_cpp(t,...)
//'     (predict(object, t+eps) - 2*object + predict(object, t-eps))/(eps^2)
//' }
//' x = seq(0,1,l=1001)
//' knots = c(1,3,5,7)/10
//' par(mfrow=c(2,2))
//' matplot(x,smf_cpp(x,knots=knots,derivs=1), type="l", lty=1)
//' matlines(x,smfD(x,knots=knots), lty=1, col=4:6)
//' abline(v=knots,lty=2)
//' matplot(x,smf_cpp(x,knots=knots,derivs=2), type="l", lty=1)
//' matlines(x,smfD2(x,knots=knots), lty=1, col=4:6)
//' abline(v=knots,lty=2)
//' matplot(x,smf_cpp(x,knots=knots,derivs=1)-smfD(x,knots=knots), type="l", lty=1)
//' abline(v=knots,lty=2)
//' matplot(x,smf_cpp(x,knots=knots,derivs=2)-smfD2(x,knots=knots), type="l", lty=1)
//' abline(v=knots,lty=2)
// [[Rcpp::export]]
Rcpp::NumericMatrix smf_cpp(vec x,
			    size_t df = 10,
			    Rcpp::Nullable<vec> knots = R_NilValue,
			    bool intercept=false,
			    size_t derivs = 0,
			    Rcpp::Nullable<arma::mat> Z = R_NilValue,
			    bool survPen_compatible = false) {
  if (intercept) {
    Rcpp::NumericMatrix out = crs_cpp(x,df,knots,derivs,true,survPen_compatible);
    out.attr("class") = "smf_cpp";
    return out;
  } else { // centering using QR
    Rcpp::NumericMatrix X0 = crs_cpp(x,df,knots,0,true,survPen_compatible);
    mat X = Rcpp::as<mat>(X0);
    mat Zi, S;
    if (Z.isNull()) {
      mat C = mat(sum(X,0).t());
      mat Q, R;
      qr(Q, R, C);
      Zi = Q(span::all, span(1,Q.n_cols-1));
    } else {
      Zi = Rcpp::as<arma::mat>(Z);
    } 
    if (derivs > 0) {
      // replacement
      X = Rcpp::as<mat>(crs_cpp(x,df,knots,derivs,true,survPen_compatible));
    }
    Rcpp::NumericMatrix out = Rcpp::wrap(X*Zi);
    out.attr("class") = Rcpp::CharacterVector{"smf_cpp","matrix"};
    out.attr("knots") = X0.attr("knots");
    out.attr("derivs") = derivs;
    out.attr("intercept") = intercept;
    out.attr("pen") = Zi.t() * Rcpp::as<mat>(X0.attr("pen")) * Zi;
    out.attr("Z") = Zi;
    return out;
  }
}
