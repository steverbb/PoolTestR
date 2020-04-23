// Generated by rstantools.  Do not edit by hand.

/*
    PoolTestR is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PoolTestR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PoolTestR.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.19.1
#include <stan/model/model_header.hpp>
namespace model_BayesianPoolScreen_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_BayesianPoolScreen");
    reader.add_event(20, 18, "end", "model_BayesianPoolScreen");
    return reader;
}
#include <stan_meta_header.hpp>
class model_BayesianPoolScreen : public prob_grad {
private:
        int N;
        std::vector<int> Result;
        std::vector<int> PoolSize;
public:
    model_BayesianPoolScreen(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_BayesianPoolScreen(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_BayesianPoolScreen_namespace::model_BayesianPoolScreen";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 1);
            current_statement_begin__ = 3;
            validate_non_negative_index("Result", "N", N);
            context__.validate_dims("data initialization", "Result", "int", context__.to_vec(N));
            Result = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("Result");
            pos__ = 0;
            size_t Result_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < Result_k_0_max__; ++k_0__) {
                Result[k_0__] = vals_i__[pos__++];
            }
            size_t Result_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < Result_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "Result[i_0__]", Result[i_0__], 0);
                check_less_or_equal(function__, "Result[i_0__]", Result[i_0__], 1);
            }
            current_statement_begin__ = 4;
            validate_non_negative_index("PoolSize", "N", N);
            context__.validate_dims("data initialization", "PoolSize", "int", context__.to_vec(N));
            PoolSize = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("PoolSize");
            pos__ = 0;
            size_t PoolSize_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < PoolSize_k_0_max__; ++k_0__) {
                PoolSize[k_0__] = vals_i__[pos__++];
            }
            size_t PoolSize_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < PoolSize_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "PoolSize[i_0__]", PoolSize[i_0__], 0);
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 7;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_BayesianPoolScreen() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 7;
        if (!(context__.contains_r("p")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable p missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("p");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "p", "double", context__.to_vec());
        double p(0);
        p = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0, 1, p);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable p: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 7;
            local_scalar_t__ p;
            (void) p;  // dummy to suppress unused var warning
            if (jacobian__)
                p = in__.scalar_lub_constrain(0, 1, lp__);
            else
                p = in__.scalar_lub_constrain(0, 1);
            // transformed parameters
            current_statement_begin__ = 10;
            validate_non_negative_index("ps", "N", N);
            std::vector<local_scalar_t__> ps(N, local_scalar_t__(0));
            stan::math::initialize(ps, DUMMY_VAR__);
            stan::math::fill(ps, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 11;
            for (int n = 1; n <= N; ++n) {
                current_statement_begin__ = 12;
                stan::model::assign(ps, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            (1 - pow((1 - p), get_base1(PoolSize, n, "PoolSize", 1))), 
                            "assigning variable ps");
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 10;
            size_t ps_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < ps_k_0_max__; ++k_0__) {
                if (stan::math::is_uninitialized(ps[k_0__])) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: ps" << "[" << k_0__ << "]";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable ps: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            size_t ps_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < ps_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "ps[i_0__]", ps[i_0__], 0);
                check_less_or_equal(function__, "ps[i_0__]", ps[i_0__], 1);
            }
            // model body
            current_statement_begin__ = 16;
            lp_accum__.add(beta_log(p, 1, 1));
            current_statement_begin__ = 17;
            lp_accum__.add(bernoulli_log(Result, ps));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("p");
        names__.push_back("ps");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_BayesianPoolScreen_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double p = in__.scalar_lub_constrain(0, 1);
        vars__.push_back(p);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 10;
            validate_non_negative_index("ps", "N", N);
            std::vector<double> ps(N, double(0));
            stan::math::initialize(ps, DUMMY_VAR__);
            stan::math::fill(ps, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 11;
            for (int n = 1; n <= N; ++n) {
                current_statement_begin__ = 12;
                stan::model::assign(ps, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            (1 - pow((1 - p), get_base1(PoolSize, n, "PoolSize", 1))), 
                            "assigning variable ps");
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 10;
            size_t ps_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < ps_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "ps[i_0__]", ps[i_0__], 0);
                check_less_or_equal(function__, "ps[i_0__]", ps[i_0__], 1);
            }
            // write transformed parameters
            if (include_tparams__) {
                size_t ps_k_0_max__ = N;
                for (size_t k_0__ = 0; k_0__ < ps_k_0_max__; ++k_0__) {
                    vars__.push_back(ps[k_0__]);
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    static std::string model_name() {
        return "model_BayesianPoolScreen";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "p";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t ps_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < ps_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "ps" << '.' << k_0__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "p";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t ps_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < ps_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "ps" << '.' << k_0__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_BayesianPoolScreen_namespace::model_BayesianPoolScreen stan_model;
#endif
