/*! \file RadMom1D_Single_Block.h
  \brief Header file defining template specializations for 1D RadMom subroutines. */

#ifndef _RADMOM_1D_SINGLE_BLOCK_INCLUDED
#define _RADMOM_1D_SINGLE_BLOCK_INCLUDED

/* Include 1D RadMom solution header file. */
#ifndef _RADMOM_1D_RECON_INCLUDED
#include "RadMom1D_Recon.h"
#endif // _RADMOM_1D_RECON_INCLUDED

/******************************************************//**
 * Routine: ICs
 *
 * Assigns initial conditions and data to the
 * solution variables.
 *
 ********************************************************/
template<class cState, class pState>
void ICs(RadMom1D_UniformMesh<cState, pState> *Soln,
         RadMom1D_Input_Parameters<cState, pState> &IP) {
  int i;
  int TC;
  double x;
  double Temperature = IP.Temperature;

  Medium1D_State MM(IP.Mo);
  cState UU(IP.Uo);

  TC = IP.Number_of_Cells_Idir+2*Soln[0].Nghost; // total number of cells

  /* Assign the initial data for the IVP of interest. */

  switch(IP.i_ICs) {
    case IC_UNIFORM :
      for ( i = 0 ; i <= TC-1 ; ++i ) {
        Soln[i].U = UU;
        Soln[i].W = Soln[i].W.U();
      } /* endfor */
      break;
    case IC_PARALLEL_PLATES :
      //----------------------------------------------------
      // Parallel Plates ICs.
      //----------------------------------------------------
      for ( i = 0 ; i <= TC-1 ; ++i ) {
        x = (Soln[i].X.x - IP.X_Min) / (IP.X_Max - IP.X_Min);
        x = max( x, 0.0 );
        x = min( x, 1.0 );

        // parabolic H2O-N2 mixture
        if (IP.Case==2) {
          MM.kappa() = 4.0*(1.0-x)*x;
        }

        if (IP.Case==1 || IP.Case==2) {
          Temperature = 300.0 + 500.0*(1.0-cos(2*PI*x));
          IP.EastWallTemp = 300.0;
          IP.WestWallTemp = 300.0;
        }

        // boundary layer temperature profile
        if (IP.Case == 3 || IP.Case == 4) {
          Temperature = 1300.0 + 350.0*cos(PI*x) - 650.0*pow(cos(PI*x),2);

          IP.EastWallTemp = 300.0;
          IP.WestWallTemp = 1000.0;
        }

        // Set initial conditions for conserved and primitive solution vectors
        Soln[i].U = UU;
        Soln[i].W = Soln[i].U.W();

        // Set gas radiative properties
        MM.SetState(Temperature);
        Soln[i].M = MM;
      }
      break;
    default:
      for ( i = 0 ; i <= TC-1 ; ++i ) {
        Soln[i].U = UU;
        Soln[i].W = Soln[i].W.U();
      } /* endfor */
      break;
  } /* endswitch */
}

/******************************************************//**
 * Routine: BCs
 *
 * Assigns boundary conditions and data (Ghost cells) to the
 * solution variables.
 *
 ********************************************************/
template<class cState, class pState>
void BCs(RadMom1D_UniformMesh<cState, pState> *Soln,
         RadMom1D_Input_Parameters<cState, pState> &IP) {

  pState Wl, Wr;
  int ICl, ICu;
  ICl = Soln[0].ICl;
  ICu = Soln[0].ICu;
  int Nghost = Soln[0].Nghost;

  switch(IP.BC_West){
    /* Apply periodic BCs for some specific problems */
    case BC_GRAY_WALL:
      Wl.Gray_Wall(Soln[ICl].W,
                   IP.WestWallTemp,
                   IP.WestWallEmiss,
                   -1);
      for (int i = 0; i < Nghost; i++) {
        Soln[(ICl - 1) - i].W = Wl;
        Soln[(ICl - 1) - i].U = Soln[(ICl - 1) - i].W.U();
      }
      break;
    case BC_CHARACTERISTIC:
      // if (IP.BCs_Enforcement == WEAK_BCS)  {
        Wl.Characteristic(Soln[ICl].W,
                          IP.WestWallTemp,
                          IP.WestWallEmiss,
                          -1);
      // } else if (IP.BCs_Enforcement == STRONG_BCS) {
      //   Wl.Gray_Wall(Soln[ICl].W,
      //                IP.WestWallTemp,
      //                IP.WestWallEmiss,
      //                -1);
      // }
      for (int i = 0; i < Nghost; i++) {
        Soln[(ICl - 1) - i].W = Wl;
        Soln[(ICl - 1) - i].U = Soln[(ICl - 1) - i].W.U();
      }
      break;
    case BC_PARTIAL_MOMENTS:
            Wl.PartialMoments_n(Soln[ICl].W,
                                IP.EastWallTemp,
                                IP.EastWallEmiss,
                                -1);
      break;

    case BC_PARTIAL_FLUX:
            // Wl.PartialFlux_n(Wr,
                             // IP.EastWallTemp,
                             // IP.EastWallEmiss);

      break;

    /* By default, constant extrapolation boundary
     *conditions are applied at either end of the mesh. */
    default:
      for (int i = 0; i < Nghost; ++i){
        // left end
        Soln[(ICl - 1) - i].U = Soln[ICl].U;
        Soln[(ICl - 1) - i].W = Soln[ICl].W;
      }
      break;
  }

  switch(IP.BC_East){
    /* Apply periodic BCs for some specific problems */
    case BC_GRAY_WALL:
      Wr.Gray_Wall(Soln[ICu].W,
                   IP.EastWallTemp,
                   IP.EastWallEmiss,
                   1);

      for (int i = 0; i < Nghost; i++) {
        Soln[(ICu + 1) + i].W = Wr;
        Soln[(ICu + 1) + i].U = Soln[(ICu + 1) + i].W.U();
      }
      break;
    case BC_CHARACTERISTIC:
      // if (IP.BCs_Enforcement == WEAK_BCS)  {
        Wr.Characteristic(Soln[ICu].W,
                          IP.EastWallTemp,
                          IP.EastWallEmiss,
                          1);
      // } else if (IP.BCs_Enforcement == STRONG_BCS) {
      //   Wr.Gray_Wall(Soln[ICu].W,
      //                IP.EastWallTemp,
      //                IP.EastWallEmiss,
      //                1);
      // }
      for (int i = 0; i < Nghost; i++) {
        Soln[(ICu + 1) + i].W = Wr;
        Soln[(ICu + 1) + i].U = Soln[(ICu + 1) + i].W.U();
      }
      break;
    case BC_PARTIAL_MOMENTS:
      Wr.PartialMoments_n(Soln[ICu].W,
                          IP.WestWallTemp,
                          IP.WestWallEmiss,
                          1);
      break;

    case BC_PARTIAL_FLUX:

      break;

    /* By default, constant extrapolation boundary
     *conditions are applied at either end of the mesh. */
    default:
      for (int i = 0; i < Nghost; ++i){
        // right end
        Soln[(ICu + 1) + i].U = Soln[ICu].U;
        Soln[(ICu + 1) + i].W = Soln[ICu].W;
      }
      break;
  }
}

/******************************************************//**
 * Routine: CFL
 *
 * Determines the allowable global and local time steps
 * (for explicit RadMom time stepping scheme) according
 * to the Courant-Friedrichs-Lewy condition.
 *
 ********************************************************/
template<class cState, class pState>
double CFL(RadMom1D_UniformMesh<cState, pState> *Soln) {

    int i;
    double dtMin;
    int NUM_VAR_RADMOM1D = Soln[0].W.NumVar();
    pState lam_x;
    double lam_max;

    /* Determine local and global time steps. */

    dtMin = MILLION;

    for ( i = Soln[0].ICl; i <= Soln[0].ICu ; ++i ) {
      // Compute wavespeed for the given cell
      lam_x = Soln[i].W.lambda_x();
      lam_max = fabs(lam_x[1]);
      for (int k = 1; k <= NUM_VAR_RADMOM1D-1; ++k){
        if (fabs(lam_max) < fabs(lam_x[k+1])){
          lam_max = fabs(lam_x[k+1]);
        }
      }

      Soln[i].dt = Soln[i].X.dx/lam_max;

      dtMin = min(dtMin, Soln[i].dt);
    } /* endfor */

    /* Return the global time step. */

    return (dtMin);

}

/********************************************************
 * Routine: L1_Norm_Residual                            *
 *                                                      *
 * Determines the L1-norm of the solution residual for  *
 * the specified quadrilateral solution block.          *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
template<class cState, class pState>
double L1_Norm_Residual(RadMom1D_UniformMesh<cState, pState> *Soln) {

    double l1_norm(ZERO);
    int norm = Soln[0].residual_variable;
    for ( int i = Soln[0].ICl; i <= Soln[0].ICu ; ++i ) {
      l1_norm += fabs(Soln[i].dUdt[norm]);
    }
    l1_norm /= Soln[0].NCi;
    return (l1_norm);
}

/********************************************************
 * Routine: L2_Norm_Residual                            *
 *                                                      *
 * Determines the L2-norm of the solution residual for  *
 * the specified quadrilateral solution block.          *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
template<class cState, class pState>
double L2_Norm_Residual( RadMom1D_UniformMesh<cState, pState> *Soln) {
    double l2_norm(ZERO);
    int norm = Soln[0].residual_variable;
    for ( int i = Soln[0].ICl; i <= Soln[0].ICu ; ++i ) {
      l2_norm += sqr(Soln[i].dUdt[norm]);
    }
    l2_norm = sqrt(l2_norm);
    l2_norm /= Soln[0].NCi;

    return (l2_norm);
}

/********************************************************
 * Routine: Max_Norm_Residual                           *
 *                                                      *
 * Determines the maximum norm of the solution residual *
 * for the specified quadrilateral solution block.      *
 * Useful for monitoring convergence of the solution    *
 * for steady state problems.                           *
 *                                                      *
 ********************************************************/
template<class cState, class pState>
double Max_Norm_Residual(RadMom1D_UniformMesh<cState, pState> *Soln) {
    double max_norm(ZERO);
    int norm = Soln[0].residual_variable;
    for ( int i = Soln[0].ICl; i <= Soln[0].ICu ; ++i ) {
      max_norm = max(max_norm, fabs(Soln[i].dUdt[norm]));
    }

    return (max_norm);
}

/****************************************************************//**
 * Reconstructs the left and right solution states at the grid
 * boundaries.  This is needed by dUdt_Multistage_Explicit() and
 * dUdt_Residual_Eval()
 ********************************************************************/
 template<class cState, class pState>
void Reconstructed_LeftandRight_States(RadMom1D_UniformMesh<cState, pState> *Soln,
                                       pState &Wl,
                                       pState &Wr,
                                       const int &i,
                                       const RadMom1D_Input_Parameters<cState, pState> &IP)
{
    // declares
    static double dX;

    int ICl, ICu;
    ICl = Soln[0].ICl;
    ICu = Soln[0].ICu;

    // cout << "i = " << i << "  " << "phi = " << Soln[i].phi << endl;

    //---------------------------------------------------------------
    // WEST BOUNDARY
    //---------------------------------------------------------------
    if (i == ICl-1 && (IP.BC_West == BC_GRAY_WALL ||
                       IP.BC_West == BC_REFLECTION ||
                       IP.BC_West == BC_MARSHAK ||
                       IP.BC_West == BC_CHARACTERISTIC ||
                       IP.BC_West == BC_PARTIAL_FLUX ||
                       IP.BC_West == BC_PARTIAL_MOMENTS)) {

        dX = -Soln[ICl].X.dx/2.0; // dX < 0 in this case
        Wr.Reconstruct(Soln[ICl].W,
                       Soln[ICl].phi,
                       Soln[ICl].dWdx,
                       dX);

        if (IP.BC_West == BC_GRAY_WALL) {
            Wl.Gray_Wall(Wr,
                         IP.WestWallTemp,
                         IP.WestWallEmiss,
                         -1);
        } else if (IP.BC_West == BC_REFLECTION) {
            // Wl.Reflect(Wr);
        } else if (IP.BC_West == BC_MARSHAK) {
            // Wl.Marshak_n(Wr,
                         // IP.WestWallTemp,
                         // IP.WestWallEmiss,
                         // -1);
        } else if (IP.BC_West == BC_CHARACTERISTIC) {
            // if (IP.BCs_Enforcement == WEAK_BCS)  {
                Wl.Characteristic(Wr,
                                  IP.WestWallTemp,
                                  IP.WestWallEmiss,
                                  -1);
            // } else if (IP.BCs_Enforcement == STRONG_BCS) {
            //     Wl.Gray_Wall(Wr,
            //                  IP.WestWallTemp,
            //                  IP.WestWallEmiss,
            //              -1);
            // } else {
            //     cout << "BCs_Enforcement type not specified!!!" << endl;
            //     exit(0);
            // }
        } else if (IP.BC_West == BC_PARTIAL_MOMENTS) {
            Wl.PartialMoments_n(Wr,
                                IP.WestWallTemp,
                                IP.WestWallEmiss,
                                -1);
        } else if (IP.BC_West == BC_PARTIAL_FLUX) {
            // Wl.PartialFlux_n(Wr,
                             // IP.WestWallTemp,
                             // IP.WestWallEmiss);
      } /* endif */

      //---------------------------------------------------------------
      // EAST BOUNDARY
      //---------------------------------------------------------------
    } else if (i == ICu && (IP.BC_East == BC_GRAY_WALL ||
                            IP.BC_East == BC_REFLECTION ||
                            IP.BC_East == BC_MARSHAK ||
                            IP.BC_East == BC_CHARACTERISTIC ||
                            IP.BC_East == BC_PARTIAL_FLUX ||
                            IP.BC_East == BC_PARTIAL_MOMENTS)) {

        dX = Soln[ICu].X.dx/2.0; // dX > 0 in this case
        Wl.Reconstruct(Soln[ICu].W,
                       Soln[ICu].phi,
                       Soln[ICu].dWdx,
                       dX);

        if (IP.BC_East == BC_GRAY_WALL) {
            Wr.Gray_Wall(Wl,
                                      IP.EastWallTemp,
                                      IP.EastWallEmiss,
                         1);
        } else if (IP.BC_East == BC_REFLECTION) {
            // Wr.Reflect(Wl);
        } else if (IP.BC_East == BC_MARSHAK) {
            // Wr.Marshak_n(Wl,
                         // IP.EastWallTemp,
                         // IP.EastWallEmiss,
                         // 1);
        } else if (IP.BC_East == BC_CHARACTERISTIC) {
            // if (IP.BCs_Enforcement == WEAK_BCS)  {
                Wr.Characteristic(Wl,
                                  IP.EastWallTemp,
                                  IP.EastWallEmiss,
                                  1);
            // } else if (IP.BCs_Enforcement == STRONG_BCS) {
            //     Wr.Gray_Wall(Wl,
            //                  IP.EastWallTemp,
            //                  IP.EastWallEmiss,
            //              1);
            // } else {
            //     cout << "BCs_Enforcement type not specified!!!" << endl;
            //     exit(0);
            // }
        } else if (IP.BC_East == BC_PARTIAL_MOMENTS) {
            Wr.PartialMoments_n(Wl,
                                IP.EastWallTemp,
                                IP.EastWallEmiss,
                                1);
        } else if (IP.BC_East == BC_PARTIAL_FLUX) {
            // Wr.PartialFlux_n(Wl,
                             // IP.EastWallTemp,
                             // IP.EastWallEmiss);
      } /* endif */
      //---------------------------------------------------------------
      // EAST face is either a normal cell or possibly a FIXED,
      // NONE or EXTRAPOLATION boundary.
      //---------------------------------------------------------------
    } else {
        dX = Soln[i].X.dx/2.0; // dX > 0 in this case
        Wl.Reconstruct(Soln[i].W,
                       Soln[i].phi,
                       Soln[i].dWdx,
                       dX);

        dX = -Soln[i+1].X.dx/2.0; // dX < 0 in this case
        Wr.Reconstruct(Soln[i+1].W,
                       Soln[i+1].phi,
                       Soln[i+1].dWdx,
                       dX);
    } // endif - X-dir

}

/******************************************************//**
 * Routine: dUdt_explicitRadMom_upwind
 *
 * This routine updates the solution using a 1st-order
 * explicit RadMom time integration and 1st-order upwind
 * spatial discretization scheme in conjunction with
 * either the Roe or HLLE flux functions.
 *
 ********************************************************/
template<class cState, class pState>
int dUdt_explicitRadMom_upwind(RadMom1D_UniformMesh<cState, pState> *Soln,
                               const RadMom1D_Input_Parameters<cState, pState> &IP,
                               double &dtMin) {
    int i;
    cState Flux;
    pState Wl, Wr;
    int ICl, ICu;
    ICl = Soln[0].ICl;
    ICu = Soln[0].ICu;

    double CFL_Number = IP.CFL_Number;
    int Flux_Function_Type = IP.i_Flux_Function;
    int Local_Time_Stepping = IP.Local_Time_Stepping;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using the first-order
       upwind scheme with a variety of flux functions. */

    LimitedLinearReconstructionOverDomain(Soln,IP);

    Soln[ICl-1].dUdt.Vacuum();
    for ( i = ICl-1 ; i <= ICu ; ++i ) {
        Soln[i+1].dUdt.Vacuum();

        // Determine left and right states
        Reconstructed_LeftandRight_States<cState, pState>( Soln, Wl, Wr, i, IP );
        // cout << "Wl = " << Wl - Soln[i].W << endl;
        // cout << "Wr = " << Wr - Soln[i+1].W << endl;
        // Wl = Soln[i].W;
        // Wr = Soln[i+1].W;
        switch(Flux_Function_Type) {
            case FLUX_FUNCTION_ROE :
                Flux = FluxRoe<cState, pState>(Wl, Wr);
                break;
            case FLUX_FUNCTION_HLLE :
                Flux = FluxHLLE<cState, pState>(Wl, Wr);
                break;
            default:
                Flux = FluxRoe<cState, pState>(Wl, Wr);
                break;
        } /* endswitch */

        Soln[i].dUdt -= Flux/Soln[i].X.dx;
        Soln[i+1].dUdt += Flux/Soln[i+1].X.dx;

        // Include source term.
        Soln[i].dUdt += Soln[i].Evaluate_Source_Term();
    } /* endfor */

    Soln[ICu+1].dUdt.Vacuum();

    // Solution Update
    for ( i = ICl ; i <= ICu ; ++i ) {
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;

        Soln[i].U += (CFL_Number*Soln[i].dt)*Soln[i].dUdt;

        if (Soln[i].U.Unphysical_Properties() ) {
            cout << "\n " << " ERROR: Negative Density and/or Energy: \n"
                 << " node = " << i << "\n U = " << Soln[i].U << "\n W = " << Soln[i].U.W() << "\n dUdt = "
                 << Soln[i].dUdt << "\n";
            return (i);
        }

        Soln[i].W = W(Soln[i].U);
    } /* endfor */

    /* Solution successfully updated. */

    return (0);
}

#endif
