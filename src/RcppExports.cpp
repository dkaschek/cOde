// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/cOde.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// wrap_cvodes
Rcpp::NumericMatrix wrap_cvodes(Rcpp::NumericVector times, Rcpp::NumericVector states_, Rcpp::NumericVector parameters_, Rcpp::NumericVector initSens_, Rcpp::DataFrame events_, Rcpp::List settings, SEXP model_, SEXP jacobian_, SEXP sens_);
static SEXP cOde_wrap_cvodes_try(SEXP timesSEXP, SEXP states_SEXP, SEXP parameters_SEXP, SEXP initSens_SEXP, SEXP events_SEXP, SEXP settingsSEXP, SEXP model_SEXP, SEXP jacobian_SEXP, SEXP sens_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type states_(states_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type parameters_(parameters_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type initSens_(initSens_SEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type events_(events_SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type settings(settingsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type jacobian_(jacobian_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type sens_(sens_SEXP);
    rcpp_result_gen = Rcpp::wrap(wrap_cvodes(times, states_, parameters_, initSens_, events_, settings, model_, jacobian_, sens_));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP cOde_wrap_cvodes(SEXP timesSEXP, SEXP states_SEXP, SEXP parameters_SEXP, SEXP initSens_SEXP, SEXP events_SEXP, SEXP settingsSEXP, SEXP model_SEXP, SEXP jacobian_SEXP, SEXP sens_SEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(cOde_wrap_cvodes_try(timesSEXP, states_SEXP, parameters_SEXP, initSens_SEXP, events_SEXP, settingsSEXP, model_SEXP, jacobian_SEXP, sens_SEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int cOde_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Rcpp::NumericMatrix(*wrap_cvodes)(Rcpp::NumericVector,Rcpp::NumericVector,Rcpp::NumericVector,Rcpp::NumericVector,Rcpp::DataFrame,Rcpp::List,SEXP,SEXP,SEXP)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP cOde_RcppExport_registerCCallable() { 
    R_RegisterCCallable("cOde", "cOde_wrap_cvodes", (DL_FUNC)cOde_wrap_cvodes_try);
    R_RegisterCCallable("cOde", "cOde_RcppExport_validate", (DL_FUNC)cOde_RcppExport_validate);
    return R_NilValue;
}
