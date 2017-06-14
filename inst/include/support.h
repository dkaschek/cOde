#include <vector>
#include <array>
#include <stdexcept>
#include <datatypes.h>
#include <nvector/nvector_serial.h>

void checkModel(int t, int neq, const std::vector<double>& initials,
                const std::vector<double>& parameters, statesRHS* model)
{
    std::vector<double> states;

    try {
        states = model(t, initials, parameters);
    } catch (std::exception& ex) {
        forward_exception_to_r(ex);
    } catch (...) {
        ::Rf_error("First call to model failed for unknown reasons.");
    }

    // Check that the number of time derivatives matches with settings
    if (states.size() != neq) {
        ::Rf_error(
            "Mismatch between number of state initials and derivatives returned from model.");
    }
}

void* createCVodes(const Rcpp::List& settings)
{
    void* cvode_mem = nullptr;

    if (Rcpp::as<std::string>(settings["method"]) == "bdf") {
        cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    } else if (Rcpp::as<std::string>(settings["method"]) == "adams") {
        cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
    } else {
        ::Rf_error("Please choose bdf or adams as method");
    }

    if (cvode_mem == nullptr)
        ::Rf_error("Could not create the CVODES solver object.");
    return cvode_mem;
}

void cvSuccess(int flag, const std::string& msg)
{
    if (flag < CV_SUCCESS) {
        throw std::runtime_error(std::string("Initializing CVodes failed:\n") + msg);
    }
}

std::vector<Event> createEventVector(const Rcpp::DataFrame& events)
{
    // Map columns of the event data frame to vectors
    // This is so ugly as data frames on C++ side can only be accessed by column
    Rcpp::NumericVector variable = events[events.offset("var")];
    Rcpp::NumericVector value    = events[events.offset("value")];
    Rcpp::NumericVector time     = events[events.offset("time")];
    Rcpp::NumericVector method   = events[events.offset("method")];

    // Fill event vector
    std::vector<Event> eventVector;
    for (int i = 0; i < events.nrows(); ++i) {
        eventVector.push_back(
            Event{static_cast<int>(variable[i]), value[i], time[i], static_cast<int>(method[i])});
    }
    std::sort(eventVector.begin(), eventVector.end());
    std::reverse(eventVector.begin(), eventVector.end());

    return eventVector;
}

void setEvent(std::vector<Event>& events, N_Vector y, N_Vector* yS, double currentTime, int neq)
{
    Event event;

    while (events.size() > 0) {
        event = events.back();
        // Pop from the back of the reverse-sorted event vector all elements
        // belonging to the current timepoint.
        if (event.time == currentTime) {
            events.pop_back();
        } else {
            break;
        }

        // cOde hands changes to states and sensitivities in one big block.
        // To figure our if a specific change affect a state or a sensitivity
        // we check the position of the variables withing the ode system. If
        // this number is larger than the neq, it must be a sensitivity.
        if (event.variable < neq) {
            // This is a state
            switch (event.method) {
            case replace: NV_Ith_S(y, event.variable) = event.value; break;
            case add: NV_Ith_S(y, event.variable) += event.value; break;
            case multiply: NV_Ith_S(y, event.variable) *= event.value; break;
            default: ::Rf_error("Failure: Unknown event type."); break;
            }
        } else {
            // This is a sensitivity
            int variable = event.variable - neq;
            int col      = static_cast<int>(variable / neq);
            int row      = variable % neq;

            switch (event.method) {
            case replace: NV_Ith_S(yS[col], row) = event.value; break;
            case add:
                ::Rf_error("Failure: Event method add is not allowed to change sensitivities.");
                break;
            case multiply: NV_Ith_S(yS[col], row) *= event.value; break;
            default: ::Rf_error("Failure: Unknown event type."); break;
            }
        }
    }
}

inline void storeStates(const N_Vector y, arma::mat& outputStates, int neq, int t)
{
    std::copy(NV_DATA_S(y), NV_DATA_S(y) + neq, outputStates.begin_col(t));
}

inline void storeSensitivities(const N_Vector* yS, arma::mat& outputSensitivities, int neq, int Ns,
                               int t)
{
    auto itOutSens = outputSensitivities.begin_col(t);
    for (auto j = 0; j < Ns; ++j) {
        auto sensData = NV_DATA_S(yS[j]);
        std::copy(sensData, sensData + neq, itOutSens);
        std::advance(itOutSens, neq);
    }
}

void storeResult(void* cvode_mem, N_Vector y, N_Vector* yS, arma::mat& outputStates,
                 arma::mat& outputSensitivities, double tretSensitivities, int t, int neq, int Ns)
{
    // Copy states to output container
    storeStates(y, outputStates, neq, t);

    // Copy sensitivities to output container, if present
    if (yS != nullptr) {
        int flag = CVodeGetSens(cvode_mem, &tretSensitivities, yS);
        if (flag < CV_SUCCESS) {
            throw std::runtime_error("Error in CVodeGetSens: Could not extract sensitivities.");
        } else {
            storeSensitivities(yS, outputSensitivities, neq, Ns, t);
        }
    }
}

void checkIntegrationStep(int flag)
{
    if (flag < CV_SUCCESS) {
        switch (flag) {
        case CV_TOO_MUCH_WORK:
            throw std::runtime_error(
                "The solver took mxstep internal steps but could not reach tout.");
            break;
        case CV_TOO_MUCH_ACC:
            throw std::runtime_error("The solver could not satisfy the accuracy demanded "
                                     "by the user for some internal step");
            break;
        case CV_ERR_FAILURE:
            throw std::runtime_error("Error test failures occured too many times during "
                                     "one internal time step or minimum step size was "
                                     "reached");
            break;
        case CV_CONV_FAILURE:
            throw std::runtime_error("Convergence test failures occurred too many times "
                                     "during one internal time step or minimum step size "
                                     "was reached.");
            break;
        case CV_LINIT_FAIL:
            throw std::runtime_error("The linear solver’s initialization function failed.");
            break;
        case CV_LSETUP_FAIL:
            throw std::runtime_error(
                "The linear solver’s setup function failed in an unrecoverable manner");
            break;
        case CV_LSOLVE_FAIL:
            throw std::runtime_error(
                "The linear solver’s solve function failed in an unrecoverable manner");
            break;
        case CV_RHSFUNC_FAIL:
            throw std::runtime_error(
                "The right hand side function failed in an unrecoverable manner");
            break;
        case CV_FIRST_RHSFUNC_ERR:
            throw std::runtime_error("The right-hand side function failed at the first call.");
            break;
        case CV_REPTD_RHSFUNC_ERR:
            throw std::runtime_error(
                "The right-hand side function had repeated recoverable errors.");
            break;
        case CV_UNREC_RHSFUNC_ERR:
            throw std::runtime_error("The right-hand side function had a recoverable "
                                     "errors but no recovery is possible.");
            break;
        case CV_BAD_T:
            throw std::runtime_error("The time t is outside the last step taken.");
            break;
        case CV_BAD_DKY: throw std::runtime_error("The output derivative vector is NULL."); break;
        case CV_TOO_CLOSE:
            throw std::runtime_error("The output and initial times are too close to each other.");
            break;
        case CV_ILL_INPUT:
            throw std::runtime_error("Input to CVode or to its solver illegal or missing.");
            break;
        default:
            throw std::runtime_error(std::string("CVodes error code: ") + std::to_string(flag));
            break;
        }
    }
}
