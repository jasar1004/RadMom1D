/*! \file RadMom1D_Recon.h
  \brief Implementation of subroutines prototyped in RadMom1D_Recon.h file. */

#ifndef _RADMOM_1D_RECON_INCLUDED
#define _RADMOM_1D_RECON_INCLUDED

/* Include 1D RadMom solution header file. */
#ifndef _RADMOM1D_IO_INCLUDED
#include "RadMom1D_IO.h"
#endif // _RADMOM1D_IO_INCLUDED

/******************************************************//**
 * Routine: Linear_Reconstruction_LeastSquares
 *
 * Peforms the reconstruction of a limited piecewise
 * linear solution state within each cell of the
 * computational mesh.  A least squares approach is
 * used in the evaluation of the unlimited solution
 * gradients.  Several slope limiters may be used.
 *
 ********************************************************/
template<class cState, class pState>
void Linear_Reconstruction_LeastSquares(RadMom1D_UniformMesh<cState, pState> *Soln,
                                        const int Limiter) {

  int i, n, n2, n_pts, index[2];
  double u0Min, u0Max, uQuad[2], phi;
  double Dx, DxDx_ave;
  pState DU, DUDx_ave;

  int ICl, ICu;
  ICl = Soln[0].ICl;
  ICu = Soln[0].ICu;

  /* Carry out the limited solution reconstruction in
   *each cell. */

  for ( i = ICl ; i <= ICu ; ++i ) {
    n_pts = 2;
    index[0] = i-1;
    index[1] = i+1;

    DUDx_ave.Vacuum();
    DxDx_ave = ZERO;

    for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
      Dx = Soln[ index[n2] ].X.x - Soln[i].X.x;
      DU = Soln[ index[n2] ].W - Soln[i].W;
      DUDx_ave += DU*Dx;
      DxDx_ave += Dx*Dx;
    } /* endfor */

    DUDx_ave = DUDx_ave/double(n_pts);
    DxDx_ave = DxDx_ave/double(n_pts);

    Soln[i].dWdx = DUDx_ave/DxDx_ave;

    for ( n = 1 ; n <= DU.NumVar() ; ++n ) {
      u0Min = Soln[i].W[n];
      u0Max = u0Min;
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
        u0Min = min(u0Min, Soln[ index[n2] ].W[n]);
        u0Max = max(u0Max, Soln[ index[n2] ].W[n]);
      } /* endfor */

      // Compute unlimited solution
      uQuad[0] = Soln[i].W[n] - HALF*Soln[i].dWdx[n]*Soln[i].X.dx;
      uQuad[1] = Soln[i].W[n] + HALF*Soln[i].dWdx[n]*Soln[i].X.dx;

      switch(Limiter) {
        case LIMITER_ZERO :
          phi = ZERO;
          break;
        case LIMITER_ONE :
          phi = ONE;
          break;
        case LIMITER_BARTH_JESPERSEN :
          phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
        case LIMITER_VENKATAKRISHNAN :
          phi = Limiter_Venkatakrishnan(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
        case LIMITER_VENKATAKRISHNAN_CORRECTED :
          phi = Limiter_Venkatakrishnan_Modified(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
        case LIMITER_VANLEER :
          phi = Limiter_VanLeer(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
        case LIMITER_VANALBADA :
          phi = Limiter_VanAlbada(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
        default:
          phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
      } /* endswitch */

      Soln[i].phi[n] = phi;
    } /* endfor */

  } /* endfor */

  // Check to make sure reconstructed solution at each of the cell's interface is realizable
  // If not, we set the slope limiters back to zero
  // pState W_recon;
  // double dX;
  // // East boundary
  // dX = Soln[i].X.dx/2.0; // dX < 0 in this case;
  // W_recon.Reconstruct(Soln[i].W, Soln[i].phi, Soln[i].dWdx, dX);
  // if (W_recon.Unphysical_Properties()) {
  //   Soln[i].phi.Vacuum();
  // }
  //
  // // West boundary
  // dX = -Soln[i].X.dx/2.0; // dX < 0 in this case;
  // W_recon.Reconstruct(Soln[i].W, Soln[i].phi, Soln[i].dWdx, dX);
  // if (W_recon.Unphysical_Properties()) {
  //   Soln[i].phi.Vacuum();
  // }

  // Gradients and limiters at ghost cells
  // Soln[0].dWdx = Soln[1].phi^Soln[1].dWdx;
  // Soln[0].phi.Ones();

  // Soln[Number_of_Cells+1].dWdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dWdx;
  // Soln[Number_of_Cells+1].phi.Ones();

}

/******************************************************//**
 * Routine: Linear_Reconstruction_GreenGauss
 *
 * Peforms the reconstruction of a limited piecewise
 * linear solution state within each cell of the
 * computational mesh.  A Green-Gauss approach is used
 * in the evaluation of the unlimited solution
 * gradients.  Several slope limiters may be used.
 *
 ********************************************************/
template<class cState, class pState>
void Linear_Reconstruction_GreenGauss(RadMom1D_UniformMesh<cState, pState> *Soln,
                                      const int Limiter) {

  int i, n;
  double u0Min, u0Max, uQuad[2], phi;

  int ICl, ICu;
  ICl = Soln[0].ICl;
  ICu = Soln[0].ICu;

  /* Carry out the limited solution reconstruction in
   *each cell. */

  for ( i = ICl ; i <= ICu ; ++i ) {
    Soln[i].dWdx = HALF*(Soln[i+1].W - Soln[i-1].W)/Soln[i].X.dx;

    for ( n = 1 ; n <= Soln[i].W.NumVar() ; ++n ) {
      u0Min = min(Soln[i-1].W[n], Soln[i].W[n]);
      u0Min = min(u0Min, Soln[i+1].W[n]);
      u0Max = max(Soln[i-1].W[n], Soln[i].W[n]);
      u0Max = max(u0Max, Soln[i+1].W[n]);

      // Compute unlimited solution
      uQuad[0] = Soln[i].W[n] - HALF*Soln[i].dWdx[n]*Soln[i].X.dx;
      uQuad[1] = Soln[i].W[n] + HALF*Soln[i].dWdx[n]*Soln[i].X.dx;

      switch(Limiter) {
        case LIMITER_ZERO :
          phi = ZERO;
          break;
        case LIMITER_ONE :
          phi = ONE;
          break;
        case LIMITER_BARTH_JESPERSEN :
          phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
        case LIMITER_VENKATAKRISHNAN :
          phi = Limiter_Venkatakrishnan(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
        case LIMITER_VENKATAKRISHNAN_CORRECTED :
          phi = Limiter_Venkatakrishnan_Modified(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
        case LIMITER_VANLEER :
          phi = Limiter_VanLeer(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
        case LIMITER_VANALBADA :
          phi = Limiter_VanAlbada(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
        default:
          phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
          break;
      } /* endswitch */

      Soln[i].phi[n] = phi;
    } /* endfor */
  } /* endfor */

  // Check to make sure reconstructed solution at each of the cell's interface is realizable
  // If not, we set the slope limiters back to zero
  pState W_recon;
  double dX;
  // East boundary
  dX = Soln[i].X.dx/2.0; // dX < 0 in this case;
  W_recon.Reconstruct(Soln[i].W, Soln[i].phi, Soln[i].dWdx, dX);
  if (W_recon.Unphysical_Properties()) {
    Soln[i].phi.Vacuum();
  }

  // West boundary
  dX = -Soln[i].X.dx/2.0; // dX < 0 in this case;
  W_recon.Reconstruct(Soln[i].W, Soln[i].phi, Soln[i].dWdx, dX);
  if (W_recon.Unphysical_Properties()) {
    Soln[i].phi.Vacuum();
  }

  // Gradients and limiters at ghost cells
  // Soln[0].dWdx = Soln[1].phi^Soln[1].dWdx;
  // Soln[0].phi.Ones();

  // Soln[Number_of_Cells+1].dWdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dWdx;
  // Soln[Number_of_Cells+1].phi.Ones();
}

/******************************************************//**
 * Routine:  LimitedLinearReconstructionOverDomain
 *
 * Peforms the reconstruction of a limited piecewise
 * linear solution state within each cell of the
 * computational mesh. The input parameters object specifies
 * all the parameters necessary to perform the reconstruction
 * (e.g. method, limiter etc.).
 *
 ********************************************************/
template<class cState, class pState>
void LimitedLinearReconstructionOverDomain(RadMom1D_UniformMesh<cState, pState> *Soln, const RadMom1D_Input_Parameters<cState, pState> &IP){

  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(Soln,
                                     IP.i_Limiter);
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(Soln,
                                       IP.i_Limiter);
    break;
  default:
    throw runtime_error("LimitedLinearReconstructionOverDomain() ERROR: Unknown reconstruction type");
    break;
  } /* endswitch */

}

#endif // _RADMOM_1D_RECON_INCLUDED
