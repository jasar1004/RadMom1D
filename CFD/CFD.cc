/*!file CFD.cc
  \brief Basic CFD Subroutines. */

/* Include the CFD header file. */

#ifndef _CFD_INCLUDED
#include "CFD.h"
#endif // _CFD_INCLUDED

/********************************************************
 * CFD -- I/O Routines.                                 *
 ********************************************************/

/********************************************************
 * Routine: Output_Progress_L2norm                      *
 *                                                      *
 * This routine writes out progress information for a   *
 * CFD calculation, including iteration level, time,    *
 * CPU time, and residual norms to the standard output  *
 * device.                                              *
 *                                                      *
 ********************************************************/
void Output_Progress_L2norm(const int Number_of_Time_Steps,
			    const double &Time,
			    const CPUTime &CPU_Time,
			    const double &Residual_L2_Norm,
			    const int First_Step,
			    const int Frequency) {

   static const int progress_character = 0;
   static const double Ratio_Residual_L2_Norm = -1.0;

   Output_Progress_L2norm(Number_of_Time_Steps,
       			  Time,
			  CPU_Time,
			  Residual_L2_Norm,
			  Ratio_Residual_L2_Norm,
			  First_Step,
			  Frequency,
			  progress_character);
}

void Output_Progress_L2norm(const int Number_of_Time_Steps,
                            const double &Time,
                            const CPUTime &CPU_Time,
                            const double &Residual_L2_Norm,
                            const double &Ratio_Residual_L2_Norm,
                            const int First_Step,
                            const int Frequency,
                            const int progress_character) {

    cout << setprecision(6);
    if (First_Step || Number_of_Time_Steps % Frequency == 0) {
      cout << endl
	   << "  n = " << Number_of_Time_Steps
	   << "   t = " << Time 
	   << "   CPU t = " << CPU_Time.min();
      cout.setf(ios::scientific);
      cout << "   L2-norm = " << Residual_L2_Norm;
      if (Ratio_Residual_L2_Norm > 0) {
         cout << "   L2-norm ratio = " << Ratio_Residual_L2_Norm;
      } /* endif */
      cout.unsetf(ios::scientific);
      cout << endl << "  ";
    } /* endif */

    switch(progress_character) {
       case 0: 
         cout << "."; 
         break;
       case 1: 
         cout << "+"; 
         break;
       case 2: 
         cout << "-"; 
         break;
       case 3: 
         cout << "o"; 
         break;
       case 4: 
         cout << ":"; 
         break;
       case 5: 
         cout << "~"; 
         break;
       default: 
         cout << "."; 
         break;
    };
    cout.flush();

}

/********************************************************
 * CFD -- Slope Limiters.                               *
 ********************************************************/

/********************************************************
 * Routine: Limiter_BarthJespersen (Slope Limiter of    *
 *          Barth and Jespersen)                        *
 *                                                      *
 * This function returns the value of the slope or      *
 * gradient limiter according to the formulation        *
 * proposed by Barth and Jespersen (1989).  Given the   *
 * minimum and maximum values of all cell centered      *
 * values used in the reconstruction and the unlimited  *
 * values of the solution at the flux quadrature        *
 * points, the routine returns the computed value of    *
 * the limiter.                                         *
 *                                                      *
 ********************************************************/
double Limiter_BarthJespersen(double *uQuad,
                              const double &u0,
                              const double &u0Min,
	      	              const double &u0Max,
			      const int nQuad) {
    int i;
    double phi, y;

    phi = ONE;

    for (i = 0; i <= nQuad-1; i++) {
       if (uQuad[i] - u0 > TOLER) {
           phi = min(phi, min(ONE,(u0Max-u0)/(uQuad[i]-u0)));
       } else if (uQuad[i] - u0 < -TOLER) {
           phi = min(phi, min(ONE,(u0Min-u0)/(uQuad[i]-u0)));
       } else {
           phi = min(phi, ONE);
       } /* endif */
    } /* endfor */

    return(phi);

}

/********************************************************
 * Routine: Limiter_Venkatakrishnan (Slope Limiter of   *
 *          Venkatakrishnan)                            *
 *                                                      *
 * This function returns the value of the slope or      *
 * gradient limiter according to the formulation        *
 * proposed by Venkatakrishnan (1993).  Given the       *
 * minimum and maximum values of all cell centered      *
 * values used in the reconstruction and the unlimited  *
 * values of the solution at the flux quadrature        *
 * points, the routine returns the computed value of    *
 * the limiter.                                         *
 *                                                      *
 ********************************************************/
double Limiter_Venkatakrishnan(double *uQuad,
                               const double &u0,
                               const double &u0Min,
	      	               const double &u0Max,
			       const int nQuad) {
    int i;
    double phi, y;

    phi = ONE;

    for (i = 0; i <= nQuad-1; i++) {
       if (uQuad[i] - u0 > TOLER) {
           y = (u0Max-u0)/(uQuad[i]-u0);
           phi = min(phi, (sqr(y)+TWO*y)/(sqr(y)+y+TWO));
       } else if (uQuad[i] - u0 < -TOLER) {
           y = (u0Min-u0)/(uQuad[i]-u0);
           phi = min(phi, (sqr(y)+TWO*y)/(sqr(y)+y+TWO));
       } else {
           phi = min(phi, ONE);
       } /* endif */
    } /* endfor */

    return(phi);

}

/********************************************************
 * Routine: Limiter_Venkatakrishnan_Modified            *
 * (Slope Limiter of Venkatakrishnan modified by        *
 *  Michalak and Ollivier-Gooch)                        *
 *                                                      *
 * This function returns the value of the slope or      *
 * gradient limiter according to the formulation        *
 * propose by Michalak and Ollivier-Gooch (2008), which *
 * is based on Venkatakrishnan's slope limiter (1993).  *
 * Given the minimum and maximum values of all cell     *
 * centered values used in the reconstruction and the   *
 * unlimited values of the solution at the flux         *
 * quadrature points, the routine returns the computed  *
 * value of the limiter.                                *
 *                                                      *
 ********************************************************/
double Limiter_Venkatakrishnan_Modified(double *uQuad,
					const double &u0,
					const double &u0Min,
					const double &u0Max,
					const int nQuad) {
    int i;
    double phi, y;

    // Parameters used by Ollivier-Gooch in the aforementioned paper
    static double yt(1.5);
    static double A((yt-2.0)/(yt*yt*yt)); // the polynomial coefficient of the cubic term
    static double B((3.0-2*yt)/(yt*yt));  // the polynomial coefficient of the quadratic term

    phi = ONE;

    for (i = 0; i <= nQuad-1; i++) {
       if (uQuad[i] - u0 > TOLER) {
           y = (u0Max-u0)/(uQuad[i]-u0);
	   if (y < yt){
	     phi = min(phi,y*y*(A*y + B) + y);
	   } else {
	     phi = min(phi,ONE);
	   }
       } else if (uQuad[i] - u0 < -TOLER) {
           y = (u0Min-u0)/(uQuad[i]-u0);
	   if (y < yt){
	     phi = min(phi,y*y*(A*y + B) + y);
	   } else {
	     phi = min(phi,ONE);
	   }
       } else {
           phi = min(phi, ONE);
       } /* endif */
    } /* endfor */

    return(phi);

}

/********************************************************
 * Routine: Limiter_VanLeer (Slope Limiter of Van Leer) *
 *                                                      *
 * This function returns the value of the slope or      *
 * gradient limiter according to the formulation        *
 * proposed by Van Leer (1978).  Given the              *
 * minimum and maximum values of all cell centered      *
 * values used in the reconstruction and the unlimited  *
 * values of the solution at the flux quadrature        *
 * points, the routine returns the computed value of    *
 * the limiter.                                         *
 *                                                      *
 ********************************************************/
double Limiter_VanLeer(double *uQuad,
                       const double &u0,
                       const double &u0Min,
	      	       const double &u0Max,
		       const int nQuad) {
    int i;
    double phi, y;

    phi = ONE;

    for (i = 0; i <= nQuad-1; i++) {
       if (uQuad[i] - u0 > TOLER) {
           y = (u0Max-u0)/(uQuad[i]-u0);
           phi = min(phi, HALF*vanleer(ONE, y));
       } else if (uQuad[i] - u0 < -TOLER) {
           y = (u0Min-u0)/(uQuad[i]-u0);
           phi = min(phi, HALF*vanleer(ONE, y));
       } else {
           phi = min(phi, ONE);
       } /* endif */
    } /* endfor */ 

    return(phi);

}

double Limiter_Minmod(double *uQuad,
		      const double &u0,
		      const double &u0Min,
		      const double &u0Max,
		      const int nQuad) {
    
    int i;
    double phi, y;

    phi = ONE;
    
    for (i = 0; i <= nQuad-1; i++) {
       if (uQuad[i] - u0 > TOLER) {
           y = (u0Max-u0)/(uQuad[i]-u0);
           phi = min(phi, HALF*minmod(ONE, y));
       } else if (uQuad[i] - u0 < -TOLER) {
           y = (u0Min-u0)/(uQuad[i]-u0);
           phi = min(phi, HALF*minmod(ONE, y));
       } else {
           phi = min(phi, ONE);
       } /* endif */
    } /* endfor */ 
    
    return(phi);

}

/********************************************************
 * Routine: Limiter_VanAlbada (Slope Limiter of         *
 *          Van Albada)                                 *
 *                                                      *
 * This function returns the value of the slope or      *
 * gradient limiter according to the formulation        *
 * proposed by Van Albada (1982).  Given the            *
 * minimum and maximum values of all cell centered      *
 * values used in the reconstruction and the unlimited  *
 * values of the solution at the flux quadrature        *
 * points, the routine returns the computed value of    *
 * the limiter.                                         *
 *                                                      *
 ********************************************************/
double Limiter_VanAlbada(double *uQuad,
                         const double &u0,
                         const double &u0Min,
	      	         const double &u0Max,
			 const int nQuad) {
    int i;
    double phi, y;

    phi = ONE;

    for (i = 0; i <= nQuad-1; i++) {
       if (uQuad[i] - u0 > TOLER) {
           y = (u0Max-u0)/(uQuad[i]-u0);
           phi = min(phi, HALF*vanalbada(ONE, y, 0.10));
       } else if (uQuad[i] - u0 < -TOLER) {
           y = (u0Min-u0)/(uQuad[i]-u0);
           phi = min(phi, HALF*vanalbada(ONE, y, 0.10));
       } else {
           phi = min(phi, ONE);
       } /* endif */
    } /* endfor */

    return(phi);

}
