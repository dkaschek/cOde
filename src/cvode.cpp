#include <string>
#include <array>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <RcppArmadillo.h>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <cvodes/cvodes_dense.h>
#include <datatypes.h>
#include <interfaces.h>
#include <support.h>

// [[Rcpp::interfaces(r, cpp)]]

//' Solve an inital value problem with cvodes.
//' 
//'
//' @description Wrapper around the solver cvodes from the Sundials suite.
//'
//' @param times Numeric vector of time points at which integration results are returned.
//'
//' @param states_ Numeric vector of inital values for states.
//'
//' @param parameters_ Numeric vector of model parameters values.
//'
//' @param initSens_ Numeric vector of inital values for sensitivities.
//'
//' @param settings List of setting passed to cvodes. For a detailed documentation of the
//'     supported setting please check the
//'     \href{http://computation.llnl.gov/projects/sundials/sundials-software}{Sundials homepage}.
//'     Supported settings are
//'     \describe{
//'     \item{\code{"jacobian"}, bool.}{
//'     For \code{"jacobian" = TRUE}, a function returning the Jacobian matrix
//'     of the system must be provided by \option{jacobian_}.}
//'
//'     \item{\code{"method"}, string, can be \code{"bdf"} or \code{"adams"}.}{
//'     The integration method used. For "bdf" \code{CVodeCreate(CV_BDF, CV_NEWTON)}
//'     is called, for "adams" \code{CVodeCreate(CV_ADAMS, CV_NEWTON)}.}
//'
//'     \item{\code{"atol"}, a scalar or a vector.}{
//'     Specifies the absolute integration tolerance. If "atol" is scalar, each
//'     state is integrated with the same absolute tolerance. If the absolute
//'     error tolerance needs to be different for each state, "atol" can be a
//'     vector holding the tolerance for each state.}
//'
//'     \item{\code{"rtol"}, scalar.}{
//'     Relative integration error tolerance.}
//'
//'     \item{\code{"which_states"}, vector.}{
//'     Return the first \code{"which_states"}. If the model has \code{N} states,
//'     \code{which_states <= N} allows to dicard all states
//'     \code{> which_states}
//'     }
//'
//'     \item{\code{"which_observed"}, vector.}{
//'     Same as \code{"which_states"}, but for observables.
//'     }
//'
//'     \item{\code{"maxsteps" = 500}, scalar.}{
//'     Maximum number of internal steps allowed to reach the next output time.
//'     While not recommended, this test can be disabled by passing
//'     \code{"maxsteps" < 0}.}
//'
//'     \item{\code{"maxord" = 12 (adams) or 5 (bdf)}, scalar.}{
//'     Maximum order of the linear multistep method. Can only be set to smaller
//'     values than default.}
//'
//'     \item{\code{"hini" = "estimated"}, scalar.}{
//'     Inital step size.}
//'
//'     \item{\code{"hmin" = 0.0}, scalar.}{
//'     Minimum absolute step size.}
//'
//'     \item{\code{"hmax" = infinity}, scalar.}{
//'     Maximum absolute step size.}
//'
//'     \item{\code{"maxerr"} = 7, scalar.}{
//'     Permitted maximum number of failed error test per step.}
//'
//'     \item{\code{"maxnonlin" = 3}, scalar.}{
//'     Permitted nonlinear solver iterations per step.}
//'
//'     \item{\code{"maxconvfail" = 10}, scalar.}{
//'     Permitted convergence failures of the nonlinear solver per step.}
//'
//'     \item{\code{"stability"} = FALSE, bool.}{
//'     Stability limit detection for the "bdf" method.}
//'
//'     \item{\code{"positive"}, bool.}{
//'     Issue an error (and abort?) in case a state becomes smaller than
//'     \option{"minimum"}.}
//'
//'     \item{\code{"minimum"}, scalar.}{
//'     Lower bound below which a state is assumed negative and reported, in
//'     case \option{\code{"positive" = TRUE}}.}
//'
//'     \item{\code{"sensitivities"} = FALSE, bool.}{
//'     Integrate sensitivities of the dynamic system.}
//'     }
//'
//' @param model_ The address of the ode model. The address is obtained as the
//'     attribute \code{address} of \code{\link[base]{getNativeSymbolInfo}}. The
//'     signature of the model function must comply to
//'
//'     \code{std::array<std::vector<double>, 2> (const double& t, const std::vector<double>&
//'     states,
//'     const std::vector<double>& parameters, const std::vector<double>& forcings)}
//'
//'     Return vector \code{std::array<std::vector<double>, 2>}
//'     \enumerate{
//'     \item
//'     First dimension holds the increments for all
//'     states.
//'     \item
//'     Second dimension holds the observed state. Not sure what these are.
//'     }
//'
//'     Argument list
//'     \code{(const double& t, const std::vector<double>& states,
//'     const std::vector<double>& parameters, const std::vector<double>& forcings)}
//'     \describe{
//'     \item{t}{
//'     Most probably the requested time point, but I am not totally sure.}
//'
//'     \item{states}{
//'     Vector of current state values.}
//'
//'     \item{parameters}{
//'     Vector of parameters values.}
//'
//'     \item{forcings}{
//'     Vector of forcings acting on the model.}
//'     }
//'
//' @param sens_ The address of the function which returns the right-hand side of the
//'     sensitivity equations. Again, this address is the attribute \code{address} obtained
//'     from the call to \code{\link[base]{getNativeSymbolInfo}}. The list of arguments is the 
//'     same as for \option{model_}.
//'     
//' @param events_ Data.frame with columns \code{var} (index of the ODE state in the vector
//'     of states.), \code{time} (time point of the event), \code{method} (one of "replace", 
//'     "add", or "multiply"), \code{value} (value associated with the event).
//'     
//' @param jacobian_ The address of the function which returns the Jacobian matrix of the
//'     model. Again, this address is the attribute \code{address} obtained
//'     from the call to \code{\link[base]{getNativeSymbolInfo}}. The function
//'     must have the signature
//'     \code{arma::mat (const double& t, const std::vector<double>& states, const
//'     std::vector<double>& parameters, const std::vector<double>& forcings)}
//'     Returned is the Jacobian matrix as an \code{arma::mat} from the
//'     Armadillo package.
//'
//'     The list of arguments is the same as for \option{model_}.
//'
//' @details This function sets up the cvodes integrator and loop over the
//'     vector \option{times} of requested time points. On success, the states
//'     of the system are returend for these time points.
//'
//'     Write something about these observations, once you got a hold of them.
//'
//' @return Matrix with nrow = (no. timepoints) and
//'     ncol = (no. states + no. observed + [(no. states)x(no. parameters)]).
//'     [(no. states)x(no. parameters)] is only returned if sensitivity equations
//'     are calculated.
//'
//'     \describe{
//'     \item{First column}{
//'     Integration time points as given in \option{times}.}
//'
//'     \item{Column 2 to no. of states + 1}{
//'     The state for the respective time point.}
//'
//'     \item{no. of states + 1 to number of states + 1 + n. of observations}{
//'     Observation for the respective time point.}
//'     }
//'
//' @author Alejandro Morales, \email{morales.s.alejandro@@gmail.com}
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix wrap_cvodes(Rcpp::NumericVector times, Rcpp::NumericVector states_,
                                Rcpp::NumericVector parameters_, Rcpp::NumericVector initSens_,
                                Rcpp::DataFrame events_, Rcpp::List settings, SEXP model_,
                                SEXP jacobian_, SEXP sens_)
{
    // Cast function pointers
    statesRHS* model = reinterpret_cast<statesRHS*>(R_ExternalPtrAddr(model_));

    auto          isJac = Rcpp::as<bool>(settings["jacobian"]);
    statesJacRHS* jacobian =
        isJac ? reinterpret_cast<statesJacRHS*>(R_ExternalPtrAddr(jacobian_)) : nullptr;

    auto              isSens = Rcpp::as<bool>(settings["sensitivities"]);
    sensitivitiesRHS* sensitivities =
        isSens ? reinterpret_cast<sensitivitiesRHS*>(R_ExternalPtrAddr(sens_)) : nullptr;

    // Convert input to standard containers
    auto stateInits(Rcpp::as<std::vector<double>>(states_));
    auto parameters(Rcpp::as<std::vector<double>>(parameters_));
    auto initSensitivities =
        isSens ? Rcpp::as<std::vector<double>>(initSens_) : std::vector<double>();

    // Test model evaluation
    const int neq = stateInits.size();
    checkModel(times[0], neq, stateInits, parameters, model);

    // Create cvodes internal data structures
    // States
    N_Vector y = N_VNew_Serial(neq);
    // Sensitivities
    const int Ns = neq + parameters.size();
    N_Vector* yS = isSens ? N_VCloneVectorArray_Serial(Ns, y) : nullptr;

    // Copy initial states into cvode state container y.
    std::copy(stateInits.cbegin(), stateInits.cend(), NV_DATA_S(y));

    // Copy initial sensitivities into cvode sensitivity container yS.
    if (isSens) {
        auto it = initSensitivities.cbegin();
        for (int i = 0; i < Ns; ++i) {
            std::copy(it, it + neq, NV_DATA_S(yS[i]));
            advance(it, neq);
        }
    }

    // Initialize output matrix
    // As armadillo is column-major, output containers are allocated such that
    // each column refers to one time point.
    const int nTimepoints = times.size();
    arma::mat outputStates(neq, nTimepoints);

    // Store initials in output matrices.
    // States
    storeStates(y, outputStates, neq, 0);
    // Sensitivities
    const int nPar  = parameters.size();
    const int nSens = neq * (neq + nPar);
    arma::mat outputSensitivities;
    if (isSens) {
        outputSensitivities.set_size(nSens, nTimepoints);
        storeSensitivities(yS, outputSensitivities, neq, Ns, 0);
    }

    // Create event vector
    auto isEvents = events_.nrows() > 0;
    auto events   = isEvents ? createEventVector(events_) : std::vector<Event>();

    // Setup first event
    if (isEvents) {
        // If the first event happens right at the start of the simulation, we
        // alter the initials in y and yS directly and _do not_ set
        // CVodeSetStopTime(). If we would set CVodeSetStopTime(), cvodes would
        // fail, as it is not allowed to integrate across a zero time span.
        if (events.back().time == times[0]) {
            setEvent(events, y, yS, times[0], neq);
        }
    }

    //////////////////////
    // Initialize CVODE //
    //////////////////////

    // Instantiate a CVODES solver object
    void*       cvode_mem = createCVodes(settings);
    UserDataIVP data_model{neq, parameters, model, jacobian, sensitivities};

    try {
        // Set error output file
        // FIXME: Errors should go somewhere. Right now, they are simply discarded.
        int flag = CVodeSetErrFile(cvode_mem, nullptr);
        cvSuccess(flag, "Error: Setting error output file.");

        // Initialize CVODES solver object
        flag = CVodeInit(cvode_mem, CVRhsFnIf, times[0], y);
        cvSuccess(flag, "Could not Initialize CVODES solver object.");

        // Set absolute and relative tolerance for integration
        flag = CVodeSStolerances(cvode_mem, settings["rtol"], settings["atol"]);
        cvSuccess(flag, "Error on setting integration tolerance.");

        // Select linear solver CVDENSE
        flag = CVDense(cvode_mem, neq);
        cvSuccess(flag, "Could not set dense linear solver.");

        // Attache user data to CVODES memory block
        flag = CVodeSetUserData(cvode_mem, &data_model);
        cvSuccess(flag, "Failure: Attach user data.");

        // Do we supply equations for the Jacobian? If so, set them.
        if (Rcpp::as<bool>(settings["jacobian"])) {
            flag = CVDlsSetDenseJacFn(cvode_mem, CVDlsDenseJacFnIf);
            cvSuccess(flag, "Failure: Setup user-supplied Jacobian function.");
        }

        // Set maximum number of steps taken by the solver to reach next output time
        flag = CVodeSetMaxNumSteps(cvode_mem, settings["maxsteps"]);
        cvSuccess(flag, "Could not set maximum number of steps.");

        // Set maximum order of the linear multistep method
        flag = CVodeSetMaxOrd(cvode_mem, settings["maxord"]);
        cvSuccess(flag, "Error: Specifying maximum order of linear multistep method. ");

        // Set initial step size
        flag = CVodeSetInitStep(cvode_mem, settings["hini"]);
        cvSuccess(flag, "Error: Setting initial step size.");

        // Set minimum step size
        flag = CVodeSetMinStep(cvode_mem, settings["hmin"]);
        cvSuccess(flag, "Error:S etting minimum step size.");

        // Set the maximum step size
        flag = CVodeSetMaxStep(cvode_mem, settings["hmax"]);
        cvSuccess(flag, "Error: Setting maximum step size.");

        // Set the maximum number of error test fails per step
        flag = CVodeSetMaxErrTestFails(cvode_mem, settings["maxerr"]);
        cvSuccess(flag, "Error: Setting error test fails.");

        // Set the maximum number of nonlinear iterations per step
        flag = CVodeSetMaxNonlinIters(cvode_mem, settings["maxnonlin"]);
        cvSuccess(flag, "Error: Setting maximum number of nonlinear solver iterations.");

        // Set the maximum number of nonlinear solver convergence failures per step
        flag = CVodeSetMaxConvFails(cvode_mem, settings["maxconvfail"]);
        cvSuccess(flag, "Error: Setting maximum number of nonlinear solver convergence failures.");

        // Should BDF stability limit detection
        flag = CVodeSetStabLimDet(cvode_mem, Rcpp::as<bool>(settings["stability"]));
        cvSuccess(flag, "Error: Setting BDF stability limit detection.");
    } catch (std::exception& ex) {
        if (y == nullptr) {
            free(y);
        } else {
            N_VDestroy_Serial(y);
        }
        if (cvode_mem == nullptr) {
            free(cvode_mem);
        } else {
            CVodeFree(&cvode_mem);
        }
        Rcpp::stop(ex.what());
    } catch (...) {
        if (y == nullptr) {
            free(y);
        } else {
            N_VDestroy_Serial(y);
        }
        if (cvode_mem == nullptr) {
            free(cvode_mem);
        } else {
            CVodeFree(&cvode_mem);
        }
        Rcpp::stop("Unknown error on initializing the CVODES ");
    }

    //////////////////////////////
    // Initialize sensitivities //
    //////////////////////////////
    if (isSens) {
        try {
            // Switch on sensitivity calculation in cvodes
            int flag = CVodeSensInit(cvode_mem, Ns, CV_SIMULTANEOUS, CVSensRhsFnIf, yS);
            cvSuccess(flag, "Error: Switch on sensitivities.");

            // FIXME: Use scalar tolerances.
            flag = CVodeSensEEtolerances(cvode_mem);
            cvSuccess(flag, "Error: Setting sensitivity tolerances.");
        } catch (std::exception& ex) {
            if (y == nullptr) {
                free(y);
            } else {
                N_VDestroy_Serial(y);
            }
            if (yS == nullptr) {
                free(yS);
            } else {
                N_VDestroyVectorArray_Serial(yS, Ns);
            }
            if (cvode_mem == nullptr) {
                free(cvode_mem);
            } else {
                CVodeFree(&cvode_mem);
            }
            Rcpp::stop(ex.what());
        } catch (...) {
            if (y == nullptr) {
                free(y);
            } else {
                N_VDestroy_Serial(y);
            }
            if (yS == nullptr) {
                free(yS);
            } else {
                N_VDestroyVectorArray_Serial(yS, Ns);
            }
            if (cvode_mem == nullptr) {
                free(cvode_mem);
            } else {
                CVodeFree(&cvode_mem);
            }
            Rcpp::stop("C++ exception (unknown reason)");
        }
    }

    ////////////
    // Events //
    ////////////
    // Setting up earliest events happening past the simulation start.
    // As all events happening right at the start are already popped of the
    // stack by the last call to setEvents() events.back().time is now the first
    // event happening past the simulation starting point.
    if (events.size() > 0) {
        int flag = CVodeSetStopTime(cvode_mem, RCONST(events.back().time));
        cvSuccess(flag, "Failure: Set event time");
    }

    ////////////////////
    // Main time loop //
    ////////////////////
    try {
        // Prepare
        double tretStates        = 0;
        double tretSensitivities = 0;

        // In each time-step, solutions are advanced by calling CVode
        for (int t = 1; t < nTimepoints; ++t) {
            int flag = CVode(cvode_mem, times[t], y, &tretStates, CV_NORMAL);
            // Error check integration step
            checkIntegrationStep(flag);
            // Store current result
            storeResult(cvode_mem, y, yS, outputStates, outputSensitivities, tretSensitivities, t,
                        neq, Ns);

            // Handle events
            if (flag == CV_TSTOP_RETURN) {
                if (isSens)
                    flag = CVodeGetSens(cvode_mem, &tretSensitivities, yS);
                setEvent(events, y, yS, tretStates, neq);
                // Reset cvode, states and sensitivities
                int flag = CVodeReInit(cvode_mem, tretStates, y);
                // The flag signals the event directly _after_ the event time point.
                // Therefore, we need to solve for the last time step again.
                // This is done by decreasing the loop counter t by one.
                t -= 1;
                cvSuccess(flag, "Failure: CVode could not be re-initialized.");
                if (isSens) {
                    flag = CVodeSensReInit(cvode_mem, CV_SIMULTANEOUS, yS);
                    cvSuccess(flag, "Failure: Sensitivities could not be re-initialized.");
                }

                // Set stop time for the next event
                if (events.size() > 0) {
                    int flag = CVodeSetStopTime(cvode_mem, RCONST(events.back().time));
                    cvSuccess(flag, "Failure: Set event time");
                }
            }
        }
    } catch (std::exception& ex) {
        if (y == nullptr) {
            free(y);
        } else {
            N_VDestroy_Serial(y);
        }
        if (yS == nullptr) {
            free(yS);
        } else {
            N_VDestroyVectorArray_Serial(yS, Ns);
        }
        if (cvode_mem == nullptr) {
            free(cvode_mem);
        } else {
            CVodeFree(&cvode_mem);
        }
        Rcpp::stop(ex.what());
    } catch (...) {
        if (y == nullptr) {
            free(y);
        } else {
            N_VDestroy_Serial(y);
        }
        if (yS == nullptr) {
            free(yS);
        } else {
            N_VDestroyVectorArray_Serial(yS, Ns);
        }
        if (cvode_mem == nullptr) {
            free(cvode_mem);
        } else {
            CVodeFree(&cvode_mem);
        }
        Rcpp::stop("C++ exception (unknown reason)");
    }

    ////////////////////////
    // Cleanup and return //
    ////////////////////////

    // Free our resources
    y == nullptr ? free(y) : N_VDestroy_Serial(y);
    yS == nullptr ? free(yS) : N_VDestroyVectorArray_Serial(yS, Ns);
    cvode_mem == nullptr ? free(cvode_mem) : CVodeFree(&cvode_mem);

    // Prepare output and return
    arma::mat outputTime(nTimepoints, 1);
    std::copy(times.begin(), times.end(), outputTime.begin_col(0));

    arma::inplace_trans(outputStates);
    arma::inplace_trans(outputSensitivities);

    arma::mat output =
        arma::join_horiz(arma::join_horiz(outputTime, outputStates), outputSensitivities);

    return Rcpp::wrap(static_cast<arma::mat>(output));
}

/////////////////////////////////
// Discard >n states on return //
/////////////////////////////////

// // Get and check number of output states and observations
// auto noutStates = Rcpp::as<long int>(settings["which_states"]);
// auto noutObserved = Rcpp::as<long int>(settings["which_observed"]);
// // This must be checked before comparisons to *.size(), as size_type is
// // unsigned.
// if(noutStates < 0 or noutObserved < 0) {
//   ::Rf_error("Negative amount of states or observables requested");
// }
//
// if(noutStates > neq or noutObserved > first_call[1].size()) {
//   ::Rf_error("More states or observations requested than provided by the model");
// }
//
// if(noutStates == 0 and noutObserved == 0) {
//   ::Rf_error("Request at least one state or observable to be returned");
// }
