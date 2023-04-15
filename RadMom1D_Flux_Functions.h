/* RadMom1D_Flux_Functions.h. */

#ifndef _RADMOM1D_FLUX_FUNCTIONS_INCLUDED
#define _RADMOM1D_FLUX_FUNCTIONS_INCLUDED

#ifndef _RADMOM1D_STATE_FIRST_ORDER_INCLUDED
#include "RadMom1DState_First_Order.h"
#endif // _RADMOM1D_STATE_FIRST_ORDER_INCLUDED

#ifndef _RADMOM1D_STATE_THIRD_ORDER_INCLUDED
#include "RadMom1DState_Third_Order.h"
#endif // _RADMOM1D_STATE_THIRD_ORDER_INCLUDED

#ifndef _CFD_INCLUDED
#include "./CFD/CFD.h"
#endif // _CFD_INCLUDED

/*********************************************************
 * Routine: FluxRoe (Roe's flux function)                *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * the "linearized" approximate Riemann solver of Roe    *
 * for the two states.  See Roe (1981).                  *
 *                                                       *
 *********************************************************/
template<class cState, class pState>
cState FluxRoe(const pState &Wl,
	      	       const pState &Wr) {

    cState Flux, Flux_L, Flux_R;
    double Correction_Roe[Wl.NumVar()];

    // Compute the numerical flux vector corresponding to the left and right states
    Flux_L = Wl.Fx();
    Flux_R = Wr.Fx();

    // Now Compute the Roe correction factor
    Wl.Compute_Correction_Factor_Roe(Correction_Roe, Wl, Wr);

    for (int i = 1; i <= Wl.NumVar(); i++) {
        Flux[i] = Flux_R[i] + Flux_L[i];
        Flux[i] -= Correction_Roe[i-1];
        Flux[i] *= HALF;
     }

     return Flux;
}

template<class cState, class pState>
cState FluxRoe(const cState &Ul,
	      	       const cState &Ur) {
   return (FluxRoe(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLE (Harten-Lax-van Leer flux function) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * the so-called Harten-Lax-van Leer approximation for   *
 * the fluxes.  See Harten, Lax, van Leer (1983).        *
 *                                                       *
 *********************************************************/
template<class cState, class pState>
cState FluxHLLE(const pState &Wl,
	      	        const pState &Wr) {

    double wavespeed_l, wavespeed_r;
    pState Wa, lambdas_l, lambdas_r, lambdas_a;
    cState Flux, dUrl;

    /* Evaluate the Roe-average primitive solution state. */

    Wa.RoeAverage(Wl, Wr);

   /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();

    // cout << lambdas_l << endl;
    // cout << lambdas_r << endl;

    // Determine max and min wavespeeds
    double lambda_a_min = ONE, lambda_l_min = ONE, lambda_r_min = ONE;
    double lambda_a_max = -ONE, lambda_l_max = -ONE, lambda_r_max = -ONE;

    for (int i = 1; i <= Wl.NumVar(); i++) {
        lambda_a_min = min(lambda_a_min, lambdas_a[i]);
        lambda_a_max = max(lambda_a_max, lambdas_a[i]);
        lambda_l_min = min(lambda_l_min, lambdas_l[i]);
        lambda_l_max = max(lambda_l_max, lambdas_l[i]);
        lambda_r_min = min(lambda_r_min, lambdas_r[i]);
        lambda_r_max = max(lambda_r_max, lambdas_r[i]);
    }

    wavespeed_l = min( min(lambda_l_min, lambda_r_min), lambda_a_min);
    wavespeed_r = max( max(lambda_r_max, lambda_l_max), lambda_a_max);

    // cout << wavespeed_l << "  " << wavespeed_r << endl;

    /* Determine the intermediate state flux. */

    if (wavespeed_l >= ZERO) {
        Flux = Wl.Fx();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.Fx();
    } else {
        Flux = (wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx());
        Flux += (wavespeed_l*wavespeed_r)*dUrl;
        Flux /= (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

template<class cState, class pState>
cState FluxHLLE(const cState &Ul,
	      	        const cState &Ur) {
   return (FluxHLLE(Ul.W(), Ur.W()));
}

#endif /* _RADMOM1D_FLUX_FUNCTIONS_INCLUDED  */
