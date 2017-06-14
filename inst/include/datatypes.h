/*
 * Typedefs of function pointers to equations formulated in C++.
 * Struct used to pass around information needed to evaluate these functions.
 *
 * Typedefs
 * - statesRHS
 *   Right hand side of the states ode, aka the model.
 *
 * - statesJacRHS
 *   Right hand side to evaluate the Jacobian of the model.
 *
 * - sensitivitiesRHS
 *   Right hand side of the sensitivities: model derivatives wrt parameteres.
 *
 *
 * Structs
 * - UserDataIVP
 *   Struct passed to cvodes containing everything for solving an IVP.
 *
 */

#ifndef __RcppSundialsDataType_h__
#define __RcppSundialsDataType_h__

#include <vector>
#include <array>
#include <RcppArmadillo.h>

typedef std::vector<double> statesRHS(const double& t, const std::vector<double>& states,
                                      const std::vector<double>& parameters);

typedef std::vector<double> statesJacRHS(const double& t, const std::vector<double>& states,
                                         const std::vector<double>& parameters);

typedef std::vector<double> sensitivitiesRHS(const double& t, const std::vector<double>& y,
                                             const std::vector<double>& yS,
                                             const std::vector<double>& p);

typedef std::array<std::vector<double>, 2> dae_in_Cpp_stl(const double& t,
                                                          const std::vector<double>& states,
                                                          const std::vector<double>& derivatives,
                                                          const std::vector<double>& parameters,
                                                          const std::vector<double>& forcings);

// Struct containing the data to solve IVP of models described in C++.
struct UserDataIVP
{
    const long int      neq;
    std::vector<double> parameters;
    statesRHS*          states;
    statesJacRHS*       jacobian;
    sensitivitiesRHS*   sensitivities;
};

enum EventMethod
{
    replace = 1,
    add,
    multiply
};

struct Event
{
    int    variable;
    double value;
    double time;
    int    method;
};

bool operator<(const Event& lhs, const Event& rhs) { return lhs.time < rhs.time; }

// Struct that contains the data to run R models
struct data_R
{
    Rcpp::NumericVector parameters;
    Rcpp::List          forcings_data;
    int                 neq;
    Rcpp::Function      model;
    Rcpp::Function      jacobian;
};

// Struct that contains the data to run C++ (with stl) DAE models
struct data_ida_Cpp_stl
{
    std::vector<double>          parameters;
    const std::vector<arma::mat> forcings_data;
    const int                    neq;
    dae_in_Cpp_stl*              model;
    statesJacRHS*                jacobian;
};

/*
 *
 Functions to interpolate the forcings when using Rcpp type of containers
 *
 */

// Interpolate the forcings when passed as a matrix
inline double interpolate(Rcpp::NumericMatrix data, double xout)
{
    Rcpp::NumericVector x = data(Rcpp::_, 0);
    Rcpp::NumericVector y = data(Rcpp::_, 1);
    // Retrieve the range of x values where xout is located
    xout = xout + sqrt(std::numeric_limits<double>::epsilon());
    // use rule = 2 for xout > x[end]
    if (xout > *(x.end() - 1)) {
        return (*(y.end() - 1));
        // use rule = 2 for xout < x[0]
    } else if (xout < x[0]) {
        return (y[0]);
        // Otherwise do linear interpolation
    } else {
        // This is supposed to be fast if the vector is ordered...
        auto   it    = std::lower_bound(x.begin(), x.end(), xout);
        int    index = (it - x.begin()) - 1;
        double inty =
            (y[index + 1] - y[index]) / (x[index + 1] - x[index]) * (xout - x[index]) + y[index];
        return (inty);
    }
}

// Interpolate the list of forcings and return the values as a NumericVector
inline Rcpp::NumericVector interpolate_list(Rcpp::List data, double t)
{
    Rcpp::NumericVector forcings(data.size());
    for (auto i = 0; i < data.size(); i++) {
        forcings[i] = interpolate(data[i], t);
    }
    return forcings;
}

/*
 *
 Functions to interpolate the forcings when usingstl and Eigen
 *
 */

// Interpolate the forcings when passed as an Armadillo matrix
inline double interpolate(const arma::mat& data, const double& xout)
{
    // use rule = 2 when xout is too high
    if (xout >= data(data.n_rows - 1, 0)) {
        return data(data.n_rows - 1, 1);
        // use rule = 2 when xout is too low
    } else if (xout <= data(0, 0)) {
        return data(0, 1);
        // Otherwise do linear interpolation
    } else {
        // This is supposed to be fast if the first column is ordered...
        auto   it    = std::lower_bound(data.begin_col(0), data.end_col(0), xout);
        int    index = (it - data.begin_col(0)) - 1;
        double inty  = (data(index + 1, 1) - data(index, 1)) /
                          (data(index + 1, 0) - data(index, 0)) * (xout - data(index, 0)) +
                      data(index, 1);
        // Use pointer arithmetics, values are stored column-by-column
        return inty;
    }
}

// Interpolate vector of forcings and return the values as a vector<double>
inline std::vector<double> interpolate_list(const std::vector<arma::mat>& data, const double& t)
{
    std::vector<double> forcings(data.size());
    for (auto i = 0; i < data.size(); i++) {
        forcings[i] = interpolate(data[i], t);
    }
    return forcings;
}

#endif
