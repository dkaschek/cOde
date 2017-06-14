/*
 * Interface between cvodes, C, and the systems of ordinary differential
 * equations which are formulated in C++.
 *
 * Interface functions:
 * - CVRhsFnIf
 *   Interface to the odes of states.
 *
 * - CVDlsDenseJacFnIf
 *   Interface to the jacobian of states.
 *
 * - CVSensRhsFnIf
 *   Interface to the odes of sensitivities.
 *
 */

#ifndef _RCPPSUNDIALSINTERFACES_H
#define _RCPPSUNDIALSINTERFACES_H

#include <vector>
#include <array>
#include <algorithm>
#include <iterator>
#include <RcppArmadillo.h>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <datatypes.h>

/** Interfacing C++ state description with cvodes.
 *
 * This function's signature matches cvodes requested signature for the
 * function which computes the ODE right-hand side for a given value of the
 * independent variable t and state vector y. This function is documented in
 * cvs_guide.pdf under the keyword CVRhsFn.
 *
 * Internally, this function calls userData->states which is the user-supplied
 * function computing the ODE right-hand side formulated in C++.
 *
 * Parameters are named according to cvs_guide.pdf, documentation is copied.
 * \param t The current value of the independent variable.
 * \param y The current value of the dependent variable vector, y(t).
 * \param ydot The output vector f(t,y).
 * \param user_data A pointer to user data, the same as the user data parameter passed to
 * CVodeSetUserData.
 */
int CVRhsFnIf(double t, N_Vector y, N_Vector ydot, void* user_data)
{
    UserDataIVP* userData = static_cast<UserDataIVP*>(user_data);

    // Transform N_Vector to std::vector
    const auto neq = userData->neq;

    // Copy current states
    std::vector<double> statesCur(neq);
    const auto          statesData = NV_DATA_S(y);
    std::copy(statesData, statesData + neq, statesCur.begin());

    // Calculate states
    std::vector<double>states = userData->states(t, statesCur, userData->parameters);

    // Copy states into result container
    std::copy(states.cbegin(), states.cend(), NV_DATA_S(ydot));

    // Indicate success
    return 0;
}

/** Interfacing C++ jacobian description with cvodes.
 *
 * This function's signature matches cvodes requested signature for the
 * function which computes the dense Jacobian for all states. This function is
 * documented in cvs_guide.pdf under the keyword CVDlsDenseJacFn.
 *
 * Internally, this function calls userData->jacobian which is the
 * user-supplied function computing the dense Jacobian for all states formulated
 * in C++.
 *
 * Parameters are named according to cvs_guide.pdf, documentation is copied.
 * \param N Problem size, number of equations.
 * \param t The current value of the independent variable.
 * \param y The current value of the dependent variable vector, namely the predicted values of y(t).
 * \param fy The current value of the vector f(t,y).
 * \param Jac The output dense Jacobian matrix (of type DlsMat).
 * \param user_data A pointer to user data, the same as the user data parameter passed to
 * CVodeSetUserData.
 * \param tmp1 N Vectors of length N which can be used as temporary storage.
 * \param tmp2 N Vectors of length N which can be used as temporary storage.
 * \param tmp3 N Vectors of length N which can be used as temporary storage.
 */
int CVDlsDenseJacFnIf(long int N, double t, N_Vector y, N_Vector fy, DlsMat Jac, void* user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    UserDataIVP* userData = static_cast<UserDataIVP*>(user_data);

    // Transform N_Vector to std::vector
    const auto neq = userData->neq;

    // Copy current states
    std::vector<double> states(neq);
    const auto          statesData = NV_DATA_S(y);
    std::copy(statesData, statesData + neq, states.begin());

    // Calculate Jacobian
    auto jacobian = userData->jacobian(t, states, userData->parameters);

    // Copy jacobian into result container
    auto itcJac = jacobian.cbegin();
    for (int i = 0; i < neq; ++i) {
        std::copy(itcJac, itcJac + neq, DENSE_COL(Jac, i));
        std::advance(itcJac, neq);
    }

    // Indicate success.
    return 0;
}

/** Interfacing C++ sensitivity description with cvodes.
 *
 * This function's signature matches cvodes requested signature for the
 * function which computes the sensitivity right-hand side for all
 * sensitivities. This function is documented in cvs_guide.pdf under the
 * keyword CVSensRhsFn.
 *
 * Internally, this function calls userData->sensitivities which is the
 * user-supplied function computing the sensitivity right-hand side for all
 * sensitivities formulated in C++.
 *
 * Parameters are named according to cvs_guide.pdf, documentation is copied.
 * \param Ns Number of parameters for which sensitivities are computed.
 * \param t The current value of the independent variable.
 * \param y The current value of the state vector, y(t).
 * \param ydot The current value of the right-hand side of the state equations.
 * \param yS The current values of the sensitivity vectors.
 * \param ySdot The output of CVSensRhsFn. On exit it must contain the sensitivity right-hand side
 * vectors.
 * \param user_data A pointer to user data, the same as the user data parameter passed to
 * CVodeSetUserData.
 * \param tmp1 N Vectors of length N which can be used as temporary storage.
 * \param tmp2 N Vectors of length N which can be used as temporary storage.
 */
int CVSensRhsFnIf(int Ns, double t, N_Vector y, N_Vector ydot, N_Vector* yS, N_Vector* ySdot,
                  void* user_data, N_Vector tmp1, N_Vector tmp2)
{
    UserDataIVP* userData = static_cast<UserDataIVP*>(user_data);

    // Transform N_Vector to std::vector
    const auto& parameters = userData->parameters;
    const auto  neq        = userData->neq;
    const auto  npar       = parameters.size();

    // Copy current states
    std::vector<double> states(neq);
    const auto          statesData = NV_DATA_S(y);
    std::copy(statesData, statesData + neq, states.begin());

    // Copy current sensitivities
    std::vector<double> sensCur(Ns * neq);
    auto                itSensCur = sensCur.begin();
    for (auto i = 0; i < Ns; ++i) {
        const auto sensitivityData = NV_DATA_S(yS[i]);
        std::copy(sensitivityData, sensitivityData + neq, itSensCur);
        std::advance(itSensCur, neq);
    }

    // Calculate sensitivities
    std::vector<double> sensitivities = userData->sensitivities(t, states, sensCur, parameters);

    // Copy sensitivities into result container
    auto itcSens = sensitivities.cbegin();
    for (int i = 0; i < neq + npar; ++i) {
        std::copy(itcSens, itcSens + neq, NV_DATA_S(ySdot[i]));
        std::advance(itcSens, neq);
    }

    // Indicate success.
    return 0;
}

#endif
