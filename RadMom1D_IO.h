/*! \file RadMom1D_IO.h
  \brief Implementation of subroutines prototyped in RadMom1D_IO.h file. */

#ifndef _RADMOM1D_IO_INCLUDED
#define _RADMOM1D_IO_INCLUDED

/* Include 1D RadMom solution header file. */
#ifndef _RADMOM1D_MESH_INCLUDED
#include "RadMom1D_Mesh.h"
#endif // _RADMOM1D_MESH_INCLUDED

/******************************************************//**
 * Routine: Output_Tecplot
 *
 * Writes the solution to specified output stream
 * suitable for plotting with TECPLOT.
 *
 ********************************************************/
template<class cState, class pState>
void Output_Tecplot(RadMom1D_UniformMesh<cState, pState> *Soln,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    ostream &out_file) {

  int ICl, ICu;

  // Set the limits of the plotted domain
  ICl = Soln[0].ICl;
  ICu = Soln[0].ICu;

  int i;

  out_file << "TITLE = \"" /*<< cfrte1D_Name()*/ << ": 1D RadMom Solution, "
	   << "Time Step/Iteration Level = "
	   << Number_of_Time_Steps
	   << ", Time = " << Time*THOUSAND << " (ms)\"" << "\n"
	   << "VARIABLES = \"x\" \\ \n"
       << "\"kp\" \\ \n"
       << "\"E\" \\ \n"
       << "\"Fx\" \\ \n"
       << "\"Sr\" \\ \n";

  for ( i = ICl ; i <= ICu ; ++i ) {
    out_file << " " << Soln[i].X.x
             << " " << Soln[i].M.kappa()
             << " " << Soln[i].U.I0()
             << " " << Soln[i].U.I1x()
             << " " << Soln[i].U.Sr(Soln[i].M) << "\n";
  } /* endfor */

  out_file << "\n";

}

#endif /* _RADMOM1D_IO_INCLUDED  */
