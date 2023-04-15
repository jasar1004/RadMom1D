/*!\file RadMom1DState_First_Order.cc
  \brief Header file defining 1D RadMom Solution State Classes. */

#ifndef _RADMOM1D_STATE_FIRST_ORDER_INCLUDED
#include "RadMom1DState_First_Order.h"
#endif // _RADMOM1D_STATE_FIRST_ORDER_INCLUDED

// /*************************************************************
//  * RadMom1D_pState_First_Order -- Create storage and assign gas constants.*
//  *************************************************************/
template <>
int RadMom1D_pState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::closure_type = MOMENT_CLOSURE_M1;
template <>
int RadMom1D_pState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::Absorption_Model = MEDIUM1D_ABSORB_GRAY;
template <>
int RadMom1D_pState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::Scattering_Func = RADIATION_SCATTER_ISO;
template <>
double RadMom1D_pState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::c = SPEED_OF_LIGHT;
template <>
double RadMom1D_pState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::a = RADIATION_CONSTANT;
template <>
double RadMom1D_pState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::C1 = PLANCK_CONSTANT;
int         RadMom1D_pState_First_Order :: NUM_VAR_RADMOM1D_FIRST_ORDER = 0;

// /*************************************************************
//  * RadMom1D_cState -- Create storage and assign gas constants.*
//  *************************************************************/
template <>
int RadMom1D_cState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::closure_type = MOMENT_CLOSURE_M1;
template <>
int RadMom1D_cState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::Absorption_Model = MEDIUM1D_ABSORB_GRAY;
template <>
int RadMom1D_cState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::Scattering_Func = RADIATION_SCATTER_ISO;
template <>
double RadMom1D_cState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::c = SPEED_OF_LIGHT;
template <>
double RadMom1D_cState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::a = RADIATION_CONSTANT;
template <>
double RadMom1D_cState<RadMom1D_cState_First_Order,
                    RadMom1D_pState_First_Order>::C1 = PLANCK_CONSTANT;
int         RadMom1D_cState_First_Order :: NUM_VAR_RADMOM1D_FIRST_ORDER = 0;

/*********************************************************
 * Routine: Rotate                                       *
 *                                                       *
 * This function returns the solution in the local       *
 * rotated frame (clockwise).                            *
 *                                                       *
 *********************************************************/
void RadMom1D_pState_First_Order :: Rotate(const double &norm_dir) {
   RadMom1D_pState_First_Order W_rotated;

   Copy_to_W(W_rotated);

   W_rotated[1] = I0();
   W_rotated[2] = norm_dir*N1x();

   Copy(W_rotated);
}

void RadMom1D_cState_First_Order :: Rotate(const double &norm_dir) {
   RadMom1D_cState_First_Order U_rotated;

   Copy_to_U(U_rotated);

   U_rotated[1] = I0();
   U_rotated[2] = norm_dir*I1x();

   Copy(U_rotated);
}

/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged (linearized)  *
 * primitive solution state given left and right        *
 * primitive solution variables.                        *
 *                                                      *
 ********************************************************/
 void RadMom1D_pState_First_Order :: RoeAverage(const RadMom1D_pState_First_Order &Wl,
                                                const RadMom1D_pState_First_Order &Wr) {
    RadMom1D_pState_First_Order Wl_Roe, Wr_Roe, Wstar_Roe;
    RadMom1D_cState_First_Order Ul, Ur, Ustar;
    
    // Determine the left and right conserved states.
    Ul = Wl.U();
    Ur = Wr.U();
    
    Wl_Roe = U_to_W_Roe(Ul);
    Wr_Roe = U_to_W_Roe(Ur);
    Wstar_Roe.AverageStates(Wl_Roe, Wr_Roe);
    
    Ustar = W_Roe_to_U(Wstar_Roe);
    
    /* Return the Roe-averged state. */
    Copy(Ustar.W());
}
 
void RadMom1D_pState_First_Order :: AverageStates(const RadMom1D_pState_First_Order &Wl,
                                                  const RadMom1D_pState_First_Order &Wr) {
    
    RadMom1D_pState_First_Order Wstar;

    Wstar[1] = HALF*(Wl.I0() + Wr.I0());
    Wstar[2] = HALF*(Wl.N1x() + Wr.N1x());
    
    /* Return the Roe-averged state. */
    Copy(Wstar);
}

/********************************************************
 * Routine: Set_ICs                                     *
 *                                                      *
 * This function returns the initial condition state    *
 * at the left, right, upper or lower boundary or the   *
 * interior domain given an incoming radiative intensity*
 *                                                      *
 ********************************************************/
 void RadMom1D_cState_First_Order :: Set_ICs(const double &Medium_Temperature) {
     
     static Medium1D_State Mwall; // a container
     
     // set wall blackbody intensity
     Mwall.setBlackBody(Medium_Temperature);
     const double Ib_wall(Mwall.Ib());

     m_values[0] = Ib_wall*FOUR*PI;
     m_values[1] = ZERO;

     // m_values[0] = Ib_wall*TWO*PI;
     // m_values[1] = Ib_wall*PI;
}

 void RadMom1D_cState_First_Order :: Set_ICs_Intensity(const double &Ib_wall) {
     m_values[0] = Ib_wall*TWO*PI;
     m_values[1] = Ib_wall*PI;
}

/********************************************************
 * Routine: Set_BCs                                     *
 *                                                      *
 * This function returns the initial boundary conditions*
 * at the left, right, upper or lower boundary or the   *
 * interior domain given an incoming radiative intensity*
 *                                                      *
 ********************************************************/
 void RadMom1D_cState_First_Order :: Set_BCs(const double *Intensity,
                                             const double norm_dir) {
    m_values[0] = Intensity[0]*TWO*PI;
    m_values[1] = -Intensity[0]*PI;

    // m_values[0] = Intensity[0]*FOUR*PI;
    // m_values[1] = ZERO;

    Rotate(norm_dir);
}

void RadMom1D_pState_First_Order :: Gray_Wall(RadMom1D_pState_First_Order W_inner,
                                              const double &wall_temperature, 
                                              const double &wall_emissivity, 
                                              const double &norm_dir) {
    
    RadMom1D_cState_First_Order U_wall;
    static Medium1D_State Mwall; // a container
    
    // set wall blackbody intensity
    Mwall.setBlackBody(wall_temperature);
    const double Ib_wall(Mwall.Ib());
    double Iw;
    double Fx_plus, E_plus;

    //------------------------------------------------
    // for a black wall
    //------------------------------------------------
    if (wall_emissivity>MICRO) {
       Iw = wall_emissivity * Ib_wall;
    } else {
       Iw = ZERO;
    }

    //------------------------------------------------
    // For grey wall.
    //------------------------------------------------
    if ( fabs(ONE - wall_emissivity)>MICRO ) {
       cout << "Double-check this implementation for M1 and P1 Gray_Wall !!!!" << endl;

       switch (RadMom1D_pState_First_Order::closure_type){
          case MOMENT_CLOSURE_P1 :
             Fx_plus = (ONE/FOUR)*W_inner.I0() + (ONE/TWO)*W_inner.I0()*W_inner.N1x();
             break;
          case MOMENT_CLOSURE_M1 :
             W_inner.Partial_Normalized_Moments_1D(E_plus, Fx_plus);
             break;
       };
       Iw = wall_emissivity * Ib_wall;
       Iw += (ONE - wall_emissivity) * Fx_plus/PI;
    }

    U_wall.Set_BCs(&Iw, norm_dir);
    Copy(U_wall.W());
}

/********************************************************
 * Routine: PartialMoments_n                               *
 *                                                      *
 * This function returns the boundary values in the     *
 * x-direction using partial fluxes given the primitive *
 * solution variables in the cell just inside the       *
 * boundary.                                            *
 *                                                      *
 ********************************************************/
void RadMom1D_pState_First_Order :: PartialMoments_n(RadMom1D_pState_First_Order W_inner,
                                                     const double &wall_temperature, 
                                                     const double &wall_emissivity,
                                                     const double &norm_dir) {
     
     double E_plus, E_minus, Fx_minus, Fx_plus;
     RadMom1D_cState_First_Order U_partial;

     W_inner.Rotate(norm_dir);
     Gray_Wall(W_inner, wall_temperature, wall_emissivity, norm_dir);
     Rotate(norm_dir);

     switch (RadMom1D_pState_First_Order::closure_type){
        case MOMENT_CLOSURE_P1 :
           E_plus = (ONE/TWO)*W_inner.I0() + (THREE/FOUR)*W_inner.I0()*W_inner.N1x();
           Fx_plus = (ONE/FOUR)*W_inner.I0() + (ONE/TWO)*W_inner.I0()*W_inner.N1x();

           E_minus = I0();
           Fx_minus = I0()*N1x();

           m_values[0] = E_plus + E_minus;
           m_values[1] = (Fx_plus + Fx_minus)/I0();
           break;
        case MOMENT_CLOSURE_M1 :
           E_minus = I0();
           Fx_minus = I0()*N1x();

           W_inner.Partial_Normalized_Moments_1D(E_plus, Fx_plus);
           m_values[0] = E_plus + E_minus;
           m_values[1] = (Fx_plus+Fx_minus)/I0();
           break;
   };
    
    Rotate(norm_dir);
}

void RadMom1D_pState_First_Order :: Partial_Normalized_Moments_1D(double &E_plus, double &Fx_plus) {
    double zeta;

    if (fabs(N1x()) < 1.0e-3) {
        E_plus = (ONE/TWO) + (THREE/FOUR)*N1x();
        Fx_plus = (ONE/FOUR) + (ONE/TWO)*N1x();
    } else {
        zeta = xi();
        E_plus = zeta*(NINE*sqr(N1x()) - SIXTEEN) - THIRTY*sqr(N1x()) + cube(N1x()) + SIX*sqr(sqr(N1x())) + THIRTY_TWO;
        E_plus /= TWO*cube(N1x());

        Fx_plus = THREE*sqr(sqr(N1x())) + zeta*cube(N1x()) - TWO*cube(N1x()) - TWO*cube(zeta) + SIX*sqr(zeta) - EIGHT;
        Fx_plus /= TWO*sqr(N1x())*(zeta - TWO);
    }

    E_plus *= I0();
    Fx_plus *= I0();
}

/********************************************************
 * Routine: Characteristic                            *
 *                                                      *
 * This function returns the boundary values using      *
 * characteristic wavespeeds in the x-direction given   *
 * the primitive solution variables at the boundary     *
 * and just inside the boundary.                        *
 *                                                      *
 ********************************************************/
 void RadMom1D_pState_First_Order :: Characteristic(RadMom1D_pState_First_Order W_inner,
                                                     const double &wall_temperature,
                                                     const double &wall_emissivity,
                                                     const double &norm_dir) {
     double lambda_val;
     RadMom1D_pState_First_Order W_charac, W_ghost;
     RadMom1D_pState_First_Order Wstar;
     RadMom1D_cState_First_Order U_bound;
     Eigenstructure_M1 Eig_M1;
     double rc_val, lc_val;

     W_inner.Rotate(norm_dir);

     Gray_Wall(W_inner, wall_temperature, wall_emissivity, norm_dir);
     Rotate(norm_dir);
     Copy_to_W(W_ghost);

     Wstar.RoeAverage(W_inner, W_ghost);

     switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
           Setup_Eigenstructure_P1(Eig_M1);
           break;
        case MOMENT_CLOSURE_M1 :
           // Precompute eigenstructure of the non-gray M1 closure
           // if (Wstar.I0() > TOLER_RADIATIVE_DENSITY_M1 &&
           //    fabs(W_inner.I0() - W_ghost.I0()) > TOLER_RADIATIVE_DENSITY_M1) {
           //    Setup_Eigenstructure_M1(Eig_M1, W_inner, W_ghost);
           // } else {
              Wstar.Setup_Eigenstructure_M1(Eig_M1);
           // }
           break;
        default:
           cout << "Closure type not specified" << endl;
           exit(0);
           break;
     }

     // Compute characteristic variables at the boundaries
     for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++) {
        lambda_val = Eig_M1.lambdas[i];

        // cout << "Eig_M1.lambdas[0] = " << Eig_M1.lambdas[0] << "  " << "Eig_M1.lambdas[1] = " << Eig_M1.lambdas[1] << endl;

        W_charac[i+1] = ZERO;
        if (lambda_val > -TOLER) {
           // Characteristic variable at the boundary is based on the incoming solution
           // which in this case corresponds to the inner solution
           for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; j++) {
              lc_val = Eig_M1.lc_vec[i][j];
              W_charac[i+1] += lc_val * W_inner.U(j+1);
           }
        } else {
           // Then characteristic variable at the boundary is based on the outgoing solution
           // which in this case corresponds to the ghost cell solution
           for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; j++) {
              lc_val = Eig_M1.lc_vec[i][j];
              W_charac[i+1] += lc_val * W_ghost.U(j+1);
         }
        }
      }

      // /////////////////////////////////////////////////////////
      //   double temp_val;
      //   for (int i = 0; i < 2; i++) {
      //       for (int j = 0; j < 2; j++) {
      //       temp_val = 0.0;
      //       for (int k = 0; k < 2; k++) {
      //           // Compute flux Jacobian based on eigen-decomposition
      //          temp_val += Eig_M1.rc_vec[i][k] * Eig_M1.lambdas[k] * Eig_M1.lc_vec[k][j];
      //               //                     cout << "rc = " << rc_vec[k][l] << "  " <<  "lc = " << lc_vec[k][l] << "  " <<  "lambdas = " << lambdas[k][l] << "  " <<  "k = " << k << "  " <<  "l = " << l << endl;
      //           }
      //
      //           // if (j == 0 )
      //              // cout << endl;
      //           // cout << temp_val << "  ";
      //
      //           if (fabs(Eig_M1.dFdU[i][j] - temp_val) > 1.0e-6) {
      //               cout << "Eigenstructure decomposition not correct: temp_val = " << temp_val << "  " <<  "dFdU[i][j] = " << Eig_M1.dFdU[i][j] << "  " <<  "i = " << i << "  " <<  "j = " << j << endl;
      //           }
      //       }
      //   }
      //   // exit(0);

      // Now compute the solution on the boundary based the vector of characteristic variables
      // at that boundary
      for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++) {
         // Initialize the conserved variable of interest
         U_bound[i+1] = ZERO;
         for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; j++) {
            rc_val = Eig_M1.rc_vec[i][j];
            U_bound[i+1] += rc_val * W_charac[j+1];
         }
       }

       m_values[0] = U_bound.W().I0();
       m_values[1] = U_bound.W().N1x();

      Rotate(norm_dir);
}
