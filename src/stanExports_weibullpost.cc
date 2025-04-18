// Generated by rstantools.  Do not edit by hand.

#include <Rcpp.h>
using namespace Rcpp ;
#include "stanExports_weibullpost.h"

RCPP_MODULE(stan_fit4weibullpost_mod) {


    class_<rstan::stan_fit<stan_model, stan::rng_t> >("rstantools_model_weibullpost")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<stan_model, stan::rng_t> ::call_sampler)
    .method("param_names", &rstan::stan_fit<stan_model, stan::rng_t> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<stan_model, stan::rng_t> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<stan_model, stan::rng_t> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<stan_model, stan::rng_t> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<stan_model, stan::rng_t> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<stan_model, stan::rng_t> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<stan_model, stan::rng_t> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<stan_model, stan::rng_t> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<stan_model, stan::rng_t> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<stan_model, stan::rng_t> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<stan_model, stan::rng_t> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<stan_model, stan::rng_t> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<stan_model, stan::rng_t> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<stan_model, stan::rng_t> ::constrained_param_names)
    .method("standalone_gqs", &rstan::stan_fit<stan_model, stan::rng_t> ::standalone_gqs)
    ;
}
