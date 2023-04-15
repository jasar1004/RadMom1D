/*!\file RadMom1DState_Third_Order.cc
  \brief Header file defining 1D RadMom Solution State Classes. */

#ifndef _RADMOM1D_STATE_THIRD_ORDER_INCLUDED
#include "RadMom1DState_Third_Order.h"
#endif // _RADMOM1D_STATE_THIRD_ORDER_INCLUDED

// /*************************************************************
//  * RadMom1D_pState_Third_Order -- Create storage and assign gas constants.*
//  *************************************************************/
template <>
int RadMom1D_pState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::closure_type = MOMENT_CLOSURE_P3;
template <>
int RadMom1D_pState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::Absorption_Model = MEDIUM1D_ABSORB_GRAY;
template <>
int RadMom1D_pState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::Scattering_Func = RADIATION_SCATTER_ISO;
template <>
double RadMom1D_pState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::c = SPEED_OF_LIGHT;
template <>
double RadMom1D_pState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::a = RADIATION_CONSTANT;
template <>
double RadMom1D_pState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::C1 = PLANCK_CONSTANT;
int         RadMom1D_pState_Third_Order :: NUM_VAR_RADMOM1D_THIRD_ORDER = 0;

// /*************************************************************
//  * RadMom1D_cState -- Create storage and assign gas constants.*
//  *************************************************************/
template <>
int RadMom1D_cState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::closure_type = MOMENT_CLOSURE_P3;
template <>
int RadMom1D_cState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::Absorption_Model = MEDIUM1D_ABSORB_GRAY;
template <>
int RadMom1D_cState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::Scattering_Func = RADIATION_SCATTER_ISO;
template <>
double RadMom1D_cState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::c = SPEED_OF_LIGHT;
template <>
double RadMom1D_cState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::a = RADIATION_CONSTANT;
template <>
double RadMom1D_cState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order>::C1 = PLANCK_CONSTANT;
int         RadMom1D_cState_Third_Order :: NUM_VAR_RADMOM1D_THIRD_ORDER = 0;

/*********************************************************
 * Routine: Rotate                                       *
 *                                                       *
 * This function returns the solution in the local       *
 * rotated frame (clockwise).                            *
 *                                                       *
 *********************************************************/
 void RadMom1D_pState_Third_Order :: Rotate(const double &norm_dir) {
     RadMom1D_pState_Third_Order W_rotated;

     W_rotated.m_values[0] = I0();

     // Rotate flux vector
     W_rotated.m_values[1] = norm_dir*N1x();
         
     // Rotate Pressure tensor
     W_rotated.m_values[2] = sqr(norm_dir)*N2xx();
         
     // Rotate third order tensor
     W_rotated.m_values[3] = cube(norm_dir)*N3xxx();

     Copy(W_rotated);
}

 void RadMom1D_cState_Third_Order :: Rotate(const double &norm_dir) {
   RadMom1D_cState_Third_Order U_rotated;
   
   U_rotated.m_values[0] = I0();
        
   // Rotate flux vector
   U_rotated.m_values[1] = norm_dir*I1x();
        
   // Rotate Pressure tensor
   U_rotated.m_values[2] = sqr(norm_dir)*I2xx();
   
   // Rotate third order tensor
   U_rotated.m_values[3] = cube(norm_dir)*I3xxx();
       
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
void RadMom1D_pState_Third_Order :: RoeAverage(const RadMom1D_pState_Third_Order &Wl,
                                               const RadMom1D_pState_Third_Order &Wr) {
    RadMom1D_cState_Third_Order Ul, Ur, Ua;
    
    /* Determine the left and right state conservative states. */
    Ul = Wl.U();
    Ur = Wr.U();
    
    /* Determine the appropriate Roe averages. */
    /* Linearized averages for now */
    Ua.m_values[0] = HALF*(Ul.I0()+Ur.I0());
    Ua.m_values[1] = HALF*(Ul.I1x()+Ur.I1x());
    Ua.m_values[2] = HALF*(Ul.I2xx()+Ur.I2xx());
    Ua.m_values[3] = HALF*(Ul.I3xxx()+Ur.I3xxx());
    
    /* Return the Roe-averged state. */
    Copy(Ua.W());
}

/********************************************************
 * Routine: Set_ICs                                     *
 *                                                      *
 * This function returns the initial condition state    *
 * at the left, right, upper or lower boundary or the   *
 * interior domain given an incoming radiative intensity*
 *                                                      *
 ********************************************************/
 void RadMom1D_cState_Third_Order :: Set_ICs(const double &Medium_Temperature) {
    static Medium1D_State Mwall; // a container
    
    // set wall blackbody intensity
    Mwall.setBlackBody(Medium_Temperature);
    const double Ib_wall(Mwall.Ib());

    m_values[0] = Ib_wall*FOUR*PI;
    m_values[1] = ZERO;
    m_values[2] = Ib_wall*FOUR*PI/THREE;
    m_values[3] = ZERO;
    
    // m_values[0] = Ib_wall*TWO*PI;
    // m_values[1] = Ib_wall*PI;
    // m_values[2] = Ib_wall*TWO*PI/THREE;
    // m_values[3] = HALF*Ib_wall*PI;
}



void RadMom1D_cState_Third_Order :: Set_ICs_Beam(const double &Medium_Temperature) {
    static Medium1D_State Mwall; // a container
    
    // set wall blackbody intensity
    Mwall.setBlackBody(Medium_Temperature);
    const double Ib_wall(Mwall.Ib());
    
    m_values[0] = Ib_wall;
    m_values[1] = Ib_wall;
    m_values[2] = Ib_wall;
    m_values[3] = Ib_wall;
}

 void RadMom1D_cState_Third_Order :: Set_ICs_Intensity(const double &Ib_wall) {
    m_values[0] = Ib_wall*FOUR*PI;
    m_values[1] = ZERO;
    m_values[2] = Ib_wall*FOUR*PI/THREE;
    m_values[3] = ZERO;
}

/********************************************************
 * Routine: Set_BCs                                     *
 *                                                      *
 * This function returns the initial boundary conditions*
 * at the left, right, upper or lower boundary or the   *
 * interior domain given an incoming radiative intensity*
 *                                                      *
 ********************************************************/
 void RadMom1D_cState_Third_Order :: Set_BCs(const double *Intensity, 
                                             const double &norm_dir) {
     // m_values[0] = Intensity[0]*FOUR*PI;
     // m_values[1] = ZERO;
     // m_values[2] = Intensity[0]*FOUR*PI/THREE;
     // m_values[3] = ZERO;

     m_values[0] = Intensity[0]*TWO*PI;
     m_values[1] = -Intensity[0]*PI;
     m_values[2] = Intensity[0]*TWO*PI/THREE;
     m_values[3] = -Intensity[0]*PI/TWO;

     Rotate(norm_dir);
}

void RadMom1D_pState_Third_Order :: Gray_Wall(RadMom1D_pState_Third_Order W_inner,
                                              const double &wall_temperature,
                                              const double &wall_emissivity,
                                              const double &norm_dir) {
    RadMom1D_cState_Third_Order U_wall;
    static Medium1D_State Mwall; // a container

     // set wall blackbody intensity
     Mwall.setBlackBody(wall_temperature);
     const double Ib_wall(Mwall.Ib());
     double Iw;
     double Fx_plus;

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
     if ( fabs(1.0-wall_emissivity)>MICRO ) {
         cout << "Double-check this implementation for P3 Gray_Wall !!!!" << endl;

         Fx_plus = (THREE/THIRTY_TWO)*W_inner.I0() + (ONE/TWO)*W_inner.U().I1x() + (FIFTEEN/THIRTY_TWO)*W_inner.U().I2xx();

         Iw = wall_emissivity * Ib_wall;
         Iw += (ONE - wall_emissivity) * Fx_plus/PI;
     }

    U_wall.Set_BCs(&Iw, norm_dir);
    Copy(U_wall.W());
}

/********************************************************
 * Routine: Reflect                                     *
 *                                                      *
 * This function returns the reflected solution state   *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
 void RadMom1D_pState_Third_Order :: Reflect(RadMom1D_pState_Third_Order W_inner,
                                             const double &norm_dir) {

     double Ibv_val, Fx_plus;
    /* Apply the frame rotation and calculate the primitive
       solution state variables in the local rotated frame
       defined by the unit normal vector. */
    W_inner.Rotate(norm_dir);
    
    Fx_plus = (THREE/THIRTY_TWO)*W_inner.I0() + (ONE/TWO)*W_inner.U().I1x() + (FIFTEEN/THIRTY_TWO)*W_inner.U().I2xx();
    Ibv_val = Fx_plus/PI;
        
    m_values[0] = Ibv_val*FOUR*PI;
    m_values[1] = ZERO;
    m_values[2] = Ibv_val*FOUR*PI/THREE;
    m_values[3] = ZERO;

    Rotate(norm_dir);
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
void RadMom1D_pState_Third_Order :: Characteristic(RadMom1D_pState_Third_Order W_inner,
                                                   const double &wall_temperature, 
                                                   const double &wall_emissivity, 
                                                   const double &norm_dir) {
     double lambda_val;
     RadMom1D_pState_Third_Order W_charac, W_ghost, W_star;
     static RadMom1D_cState_Third_Order rc_vec, lc_vec;
     RadMom1D_cState_Third_Order U_bound;
     RadMom1D_cState_Third_Order U_inner, U_ghost;
     RadMom1D_pState_Third_Order lambda_W_star;
     Eigenstructure_P3 Eig_P3;
     double rc_val, lc_val;
     
     W_inner.Rotate(norm_dir);
     Gray_Wall(W_inner, wall_temperature, wall_emissivity, norm_dir);
     Rotate(norm_dir);
     
     Copy_to_W(W_ghost);
     
     // Compute W_star for the linearization of the flux Jacobian (based on Roe Average)
     // Note that, in the rotated frame, the inner solution corresponds to the left solution
     // and the ghost cell solution corresponds to the right solution
     W_star.RoeAverage(W_inner, W_ghost);
     
     Setup_Eigenstructure_P3(Eig_P3);

     U_inner = W_inner.U();
     U_ghost = W_ghost.U();
     lambda_W_star = W_star.lambda_x();

     // Compute characteristic variables at the boundaries
     for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++) {
        lambda_val = lambda_W_star[i+1];
             
        W_charac.m_values[i] = ZERO;
        if (lambda_val > TOLER) {
            // Characteristic variable at the boundary is based on the incoming solution
            // which in this case corresponds to the inner solution
            for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; j++) {
                lc_val = Eig_P3.lc_vec[i][j];
                W_charac.m_values[i] += lc_val*U_inner[j+1];
            }
        } else if (lambda_val <= TOLER) {
            // Then characteristic variable at the boundary is based on the outgoing solution
            // which in this case corresponds to the ghost cell solution
            for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; j++) {
                lc_val = Eig_P3.lc_vec[i][j];
                W_charac.m_values[i] += lc_val*U_ghost[j+1];
            }
        }

        // Now compute the solution on the boundary based the vector of characteristic variables
        // at that boundary
        for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++) {
            // Initialize the primitie variable of interest
            U_bound.m_values[i] = ZERO;
            for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; j++) {
                rc_val = Eig_P3.rc_vec[i][j];
                U_bound.m_values[i] += rc_val * W_charac.m_values[j];
            }
        }
        
        m_values[0] = U_bound.W().I0();
        m_values[1] = U_bound.W().N1x();
        m_values[2] = U_bound.W().N2xx();
        m_values[3] = U_bound.W().N3xxx();
     }
     
     Rotate(norm_dir);
}

/********************************************************
 * Routine: PartialFlux_n                               *
 *                                                      *
 * This function returns the boundary values in the     *
 * x-direction using partial fluxes given the primitive *
 * solution variables in the cell just inside the       *
 * boundary.                                            *
 *                                                      *
 ********************************************************/
void RadMom1D_pState_Third_Order :: PartialMoments_n(RadMom1D_pState_Third_Order W_inner,
                                                     const double &wall_temperature, 
                                                     const double &wall_emissivity, 
                                                     const double &norm_dir) {
     double E_plus, E_minus, Fx_minus, Fx_plus;
     double Pxx_plus, Pxx_minus;
     double Qxxx_plus, Qxxx_minus;
     
     W_inner.Rotate(norm_dir);
     Gray_Wall(W_inner, wall_temperature, wall_emissivity, norm_dir);
     Rotate(norm_dir);
     
     E_minus = I0();
     Fx_minus = I0()*N1x();
     Pxx_minus = I0()*N2xx();
     Qxxx_minus = I0()*N3xxx();
                 
     E_plus = (ONE/TWO)*W_inner.I0() + (FORTY_FIVE/THIRTY_TWO)*W_inner.U().I1x() - (THIRTY_FIVE/THIRTY_TWO)*W_inner.U().I3xxx();
     Fx_plus = (THREE/THIRTY_TWO)*W_inner.I0() + (ONE/TWO)*W_inner.U().I1x() + (FIFTEEN/THIRTY_TWO)*W_inner.U().I2xx();
     Pxx_plus = (FIVE/THIRTY_TWO)*W_inner.U().I1x() + (ONE/TWO)*W_inner.U().I2xx() + (THIRTY_FIVE/NINETY_SIX)*W_inner.U().I3xxx();
     Qxxx_plus = -(ONE/THIRTY_TWO)*W_inner.I0() + (FIFTEEN/THIRTY_TWO)*W_inner.U().I2xx()+ (ONE/TWO)*W_inner.U().I3xxx();

     m_values[0] = E_plus + E_minus;
     m_values[1] = (Fx_plus+Fx_minus)/I0();
     m_values[2] = (Pxx_plus+Pxx_minus)/I0();
     m_values[3] = (Qxxx_plus+Qxxx_minus)/I0();
    
     Rotate(norm_dir);
}

/********************************************************
 * Routine: Markshak_n                                  *
 *                                                      *
 * This function returns the boundary values using      *
 * standard Markshak condition from P1 method, given    *
 * the primitive solution variables at the boundary     *
 * and just inside the boundary.                        *
 *                                                      *
 ********************************************************/
void RadMom1D_pState_Third_Order :: Marshak_n(RadMom1D_pState_Third_Order W_inner,
                                               const double &wall_temperature, 
                                              const double &wall_emissivity, 
                                              const double &norm_dir) {
    double Fx_plus, Fx_minus;
    double Qxxx_plus, Qxxx_minus;
    
    W_inner.Rotate(norm_dir);
    Gray_Wall(W_inner, wall_temperature, wall_emissivity, norm_dir);
    Rotate(norm_dir);
    
    Fx_minus = I0()*N1x();
    Qxxx_minus = I0()*N3xxx();
        
    Fx_plus = (THREE/THIRTY_TWO)*W_inner.I0() + (ONE/TWO)*W_inner.U().I1x() + (FIFTEEN/THIRTY_TWO)*W_inner.U().I2xx();
    Qxxx_plus = -(ONE/THIRTY_TWO)*W_inner.I0() + (FIFTEEN/THIRTY_TWO)*W_inner.U().I2xx()+ (ONE/TWO)*W_inner.U().I3xxx();
       
    m_values[0] = FOUR*((Fx_plus - Fx_minus) - (Qxxx_plus - Qxxx_minus));
    m_values[1] = (Fx_plus + Fx_minus)/I0();
    m_values[2] = (FOUR/(FIFTEEN*I0()))*((Fx_plus - Fx_minus) + THREE*(Qxxx_plus - Qxxx_minus));
    m_values[3] = (Qxxx_plus + Qxxx_minus)/I0();
    
    Rotate(norm_dir);
}

/********************************************************
 * Routine: PartialFlux_n                               *
 *                                                      *
 * This function returns the boundary values in the     *
 * x-direction using partial fluxes given the primitive *
 * solution variables in the cell just inside the       *
 * boundary.                                            *
 *                                                      *
 ********************************************************/
void RadMom1D_pState_Third_Order :: PartialFlux_n (RadMom1D_pState_Third_Order W_inner,
                                                   const double &wall_temperature, 
                                                   const double &wall_emissivity, 
                                                   const double &norm_dir) {
    double Fx_plus, Fx_minus, Pxx_plus, Pxx_minus;
    double Qxxx_plus, Qxxx_minus;
    
    W_inner.Rotate(norm_dir);
    Gray_Wall(W_inner, wall_temperature, wall_emissivity, norm_dir);
    Rotate(norm_dir);
    
    Fx_minus = I0()*N1x();
    Pxx_minus = I0()*N2xx();
    Qxxx_minus = I0()*N3xxx();
        
    Fx_plus = (THREE/THIRTY_TWO)*W_inner.I0() + (ONE/TWO)*W_inner.U().I1x() + (FIFTEEN/THIRTY_TWO)*W_inner.U().I2xx();
    Pxx_plus = (FIVE/THIRTY_TWO)*W_inner.U().I1x() + (ONE/TWO)*W_inner.U().I2xx() + (THIRTY_FIVE/NINETY_SIX)*W_inner.U().I3xxx();
    Qxxx_plus = -(ONE/THIRTY_TWO)*W_inner.I0() + (FIFTEEN/THIRTY_TWO)*W_inner.U().I2xx()+ (ONE/TWO)*W_inner.U().I3xxx();
        
    m_values[0] = FIFTEEN*(Pxx_plus + Pxx_minus) - SIXTEEN*(Qxxx_plus - Qxxx_minus);
    m_values[1] = (Fx_plus + Fx_minus)/I0();
    m_values[2] = (Pxx_plus + Pxx_minus)/I0();
    m_values[3] = (Qxxx_plus + Qxxx_minus)/I0();
        
    Rotate(norm_dir);
    
}
