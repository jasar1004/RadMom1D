/*!\file RadMom1D_Mesh.h
  \brief Header file defining 1D RadMom Solution Classes. */

#ifndef _RADMOM1D_MESH_INCLUDED
#define _RADMOM1D_MESH_INCLUDED

/* Include 1D RadMom state and cell header files. */

#ifndef _RADMOM1D_STATE_INCLUDED
#include "RadMom1D_Flux_Functions.h"
#endif // _RADMOM1D_STATE_INCLUDED

#ifndef _CELL1D_INCLUDED
#include "./Grid/Cell1D.h"
#endif // _CELL1D_INCLUDED

#ifndef _RADMOM1D_INPUT_INCLUDED
#include "RadMom1DInput.h"
#endif // _RADMOM1D_INPUT_INCLUDED

template<class cState, class pState>
class RadMom1D_UniformMesh;

/* Relational operators. */
template<class cState, class pState>
int operator ==(const RadMom1D_UniformMesh<cState, pState> &Soln1,
			 const RadMom1D_UniformMesh<cState, pState> &Soln2);

template<class cState, class pState>
int operator !=(const RadMom1D_UniformMesh<cState, pState> &Soln1,
			 const RadMom1D_UniformMesh<cState, pState> &Soln2);

/* Input-output operators. */
template<class cState, class pState>
ostream &operator << (ostream &out_file,
			       const RadMom1D_UniformMesh<cState, pState> &Soln);

template<class cState, class pState>
istream &operator >> (istream &in_file,
			       RadMom1D_UniformMesh<cState, pState> &Soln);

/* Define the classes. */

/******************************************************//**
 * Class: RadMom1D_UniformMesh
 * 
 * @brief Class definition of the 1D RadMom computational cell (i.e. solution + geometry)
 *                                                      
 * Member functions                                     
 * -     W       -- Return primitive solution state.     
 * -     U       -- Return conserved solution state.     
 * -     X       -- Return cell geometry.                
 * -     dt      -- Return local time step.              
 * -   dUdt      -- Return the solution residual.        
 * -   dUdx      -- Return the unlimited solution        
 * -                gradient.                            
 * -    phi      -- Return the solution slope limiters.  
 * -     Uo      -- Return initial solution state.       
 * - Nghost      -- Number of ghost cells (!= 1 for      
 * -                high-order)
 *                                                      
 * Member operators \n                                    
 *      S -- a 1D RadMom solution
 *                                                      
 * - S = S;                                               
 * - S = S + S;                                           
 * - S = S - S;                                           
 * - S = +S;                                              
 * - S = -S;                                              
 * - S += S;                                              
 * - S -= S;                                              
 * - S == S;                                              
 * - S != S;                                              
 * - cout << S; (output function)                         
 * - cin  >> S; (input function)                          
 *                                                      
 ********************************************************/

template<class cState, class pState>
class RadMom1D_UniformMesh{
public:
  // typedef HighOrder1D<pState> HighOrderType; //!< high-order variable data type
  typedef pState SOLN_pSTATE;                //!< primitive solution state type
  typedef cState SOLN_cSTATE;                //!< conserved solution state type

  pState    W;   //!< Primitive solution state.
  cState    U;   //!< Conserved solution state.
  Medium1D_State   M; ; //!< Participating medium state.
  Cell1D_Uniform    X;   //!< Cell geometry.
  double           dt;   //!< Local time step.
  cState dUdt;   //!< Solution residual.
  cState TotaldUdt; //!< Solution residual. Used for RK4
  pState dWdx;   //!< Unlimited solution gradient.
  pState  phi;   //!< Solution slope limiter.
  cState   Uo;   //!< Initial solution state.
  static int residual_variable; //!< Indicates which variable is used
  // Made public so can access them.

  int Nghost;            //!< Number of ghost cells(!= 1 for high-order)
  int NCi;            //!< Number of cells
  int ICl, ICu;	   //!< Indexes (Start & End)

  /* Creation, copy, and assignment constructors. */
  RadMom1D_UniformMesh(void);
  RadMom1D_UniformMesh(const RadMom1D_UniformMesh &Soln);
  RadMom1D_UniformMesh(const pState &W0, const cState &U0, const Medium1D_State &M0, const Cell1D_Uniform &X0);
  RadMom1D_UniformMesh(const pState &W0, const Medium1D_State &M0, const Cell1D_Uniform &X0);
  RadMom1D_UniformMesh(const cState &U0, const Medium1D_State &M0, const Cell1D_Uniform &X0);

  /* Evaluate source term */
  cState Evaluate_Source_Term( void );

  /* Field access */
  const double & CellCenter(void) const {return X.x;} //!< return cell center
  const double & CellDelta (void) {return X.dx;} //!< return cell delta

  /* Operating functions */
  double SolutionAtCoordinates_PWL (const double & X_Coord, const unsigned parameter) ;
    
  /* Relational operators. */
  friend int operator == <cState, pState> (const RadMom1D_UniformMesh<cState, pState> &Soln1,
			 const RadMom1D_UniformMesh<cState, pState> &Soln2);
  friend int operator != <cState, pState> (const RadMom1D_UniformMesh<cState, pState> &Soln1,
			 const RadMom1D_UniformMesh<cState, pState> &Soln2);
    
  /* Input-output operators. */
  friend ostream &operator << <cState, pState> (ostream &out_file,
			       const RadMom1D_UniformMesh<cState, pState> &Soln);
  friend istream &operator >> <cState, pState> (istream &in_file,
			       RadMom1D_UniformMesh<cState, pState> &Soln);

private:
  
};

/**********************************************************************
 * RadMom1D_UniformMesh -- Create storage for the static variables.     *
 **********************************************************************/
// Initialize residual_variable
template<class cState, class pState>
int RadMom1D_UniformMesh<cState, pState>::residual_variable = 1;

/*******************************************************
 *  RadMom1D_UniformMesh<cState, pState>::MemberFunctions()             *
 ******************************************************/
/* Constructor */
template<class cState, class pState>
inline RadMom1D_UniformMesh<cState, pState>::RadMom1D_UniformMesh(void) {
  W.Vacuum(), M.Nullify(), U.Vacuum();
  X = Cell1D_Uniform_ONE; dt = ZERO;
  dUdt.Vacuum();
  dWdx.Vacuum(); phi.Vacuum();
  Uo.Vacuum();
  residual_variable = 1;
}

template<class cState, class pState>
inline RadMom1D_UniformMesh<cState, pState>::RadMom1D_UniformMesh(const RadMom1D_UniformMesh &Soln) {
  W = Soln.W; U = Soln.U; M = Soln.M; X = Soln.X;
  dt = Soln.dt; dUdt = Soln.dUdt; dWdx = Soln.dWdx;
  phi = Soln.phi; Uo = Soln.Uo;
  NCi = Soln.NCi; Nghost = Soln.Nghost; ICl = Soln.ICl; ICu = Soln.ICu;
  residual_variable = Soln.residual_variable;
}

template<class cState, class pState>
inline RadMom1D_UniformMesh<cState, pState>::RadMom1D_UniformMesh(const pState &W0,
						const cState &U0,
						const Medium1D_State &M0,
						const Cell1D_Uniform &X0) {
  W = W0; U = U0; M = M0; X = X0;
  dt = ZERO; dUdt.Vacuum();
  dWdx.Vacuum(); phi.Vacuum();
  Uo.Vacuum();
  residual_variable = 1;
}

template<class cState, class pState>
inline RadMom1D_UniformMesh<cState, pState>::RadMom1D_UniformMesh(const pState &W0,
                                                                  const Medium1D_State &M0,
                                                                  const Cell1D_Uniform &X0) {
  W = W0; U = W0.U(); M = M0; X = X0;
  dt = ZERO; dUdt.Vacuum();
  dWdx.Vacuum(); phi.Vacuum();
  Uo.Vacuum();
  residual_variable = 1;
}

template<class cState, class pState>
inline RadMom1D_UniformMesh<cState, pState>::RadMom1D_UniformMesh(const cState &U0,
                                                                  const Medium1D_State &M0,
                                                                  const Cell1D_Uniform &X0) {
  W = U0.W(); U = U0; M = M0; X = X0;
  dt = ZERO; dUdt.Vacuum();
  dWdx.Vacuum(); phi.Vacuum();
  Uo.Vacuum();
  residual_variable = 1;
}

//! Return the solution of the piecewise limited linear reconstruction at the coordinate X_Coord,
//  for the required parameter.
template<class cState, class pState>
inline double RadMom1D_UniformMesh<cState, pState>::SolutionAtCoordinates_PWL (const double & X_Coord, const unsigned parameter) {
  SOLN_pSTATE W_Temp;
  double Distance(X_Coord - X.x);
  W_Temp = W + (phi^dWdx)*Distance;
  return W_Temp[parameter];
}

/********************************************************
 * RadMom1D_UniformMesh -- Relational operators.         *
 ********************************************************/
template<class cState, class pState>
inline int operator ==(const RadMom1D_UniformMesh<cState, pState> &Soln1,
		       const RadMom1D_UniformMesh<cState, pState> &Soln2) {
  return (Soln1.W == Soln2.W && Soln1.U == Soln2.U && Soln1.M == Soln2.M &&
	  Soln1.X == Soln2.X && Soln1.dt == Soln2.dt &&
	  Soln1.dUdt == Soln2.dUdt && Soln1.dWdx == Soln2.dWdx &&
	  Soln1.phi == Soln2.phi && Soln1.Uo == Soln2.Uo );
}

template<class cState, class pState>
inline int operator !=(const RadMom1D_UniformMesh<cState, pState> &Soln1,
		       const RadMom1D_UniformMesh<cState, pState> &Soln2) {
  return (Soln1.W != Soln2.W || Soln1.U != Soln2.U || Soln1.M != Soln2.M ||
	  Soln1.X != Soln2.X || Soln1.dt != Soln2.dt ||
	  Soln1.dUdt != Soln2.dUdt || Soln1.dWdx != Soln2.dWdx ||
	  Soln1.phi != Soln2.phi || Soln1.Uo != Soln2.Uo );
}

/********************************************************
 * RadMom1D_UniformMesh -- Input-output operators.       *
 ********************************************************/
template<class cState, class pState>
inline ostream &operator << (ostream &out_file,
			     const RadMom1D_UniformMesh<cState, pState> &Soln) {
  out_file << Soln.X << Soln.M << Soln.U;
  out_file.setf(ios::scientific);
  return (out_file);
}

template<class cState, class pState>
inline istream &operator >> (istream &in_file,
			     RadMom1D_UniformMesh<cState, pState> &Soln) {
  in_file >> Soln.X >> Soln.M >> Soln.U;
  Soln.W = Soln.U.W();
  return (in_file);
}

/**************************************************************************
 * RadMom1D_UniformMesh::Evaluate_Source_Term.                            *
 **************************************************************************/
template<class cState, class pState>
inline cState RadMom1D_UniformMesh<cState, pState>::Evaluate_Source_Term(void) {
    cState Source;

    // Include general source term.
    Source = U.S(M);

    return Source;
}

/******************************************************//**
 * Routine: Allocate
 *
 * Allocate memory for 1D RadMom equation solution.
 *
 ********************************************************/
template<class cState, class pState>
RadMom1D_UniformMesh<cState, pState>* Allocate(RadMom1D_UniformMesh<cState, pState> *Soln_ptr,
                               const RadMom1D_Input_Parameters<cState, pState> &IP) {

  int NC;                       // number of cells in the computational domain
  int Nghost; 			// number of ghost cells
  int NCi; 			// number of ghost cells

  /* Calculate the total number of computational cells */
  Nghost = IP.Number_of_Ghost_Cells;
  NCi = IP.Number_of_Cells_Idir;
  NC = IP.Number_of_Cells_Idir + 2 * Nghost;

  /* Allocate memory. */
  Soln_ptr = new RadMom1D_UniformMesh<cState, pState>[NC];

  /* Set preliminary mesh parameters */
  for (int i=0; i<= NC-1; ++i){

    // store domain indexes in each cell
    Soln_ptr[i].Nghost = Nghost;
    Soln_ptr[i].NCi = NCi;
    Soln_ptr[i].ICl = Nghost;
    Soln_ptr[i].ICu = NC - 1 - Nghost;
  }//endfor

  /* Return memory location. */

  return(Soln_ptr);
}

/******************************************************//**
 * Routine: Deallocate
 *
 * Deallocate memory for 1D RadMom equation solution.
 *
 ********************************************************/
template<class cState, class pState>
RadMom1D_UniformMesh<cState, pState>* Deallocate(RadMom1D_UniformMesh<cState, pState> *Soln_ptr) {

  /* Deallocate memory. */
  delete []Soln_ptr;
  Soln_ptr = NULL;

  /* Return memory location. */

  return(Soln_ptr);
}

/******************************************************//**
 * Routine: Grid
 *
 * Generates a uniform mesh and assign the locations of
 * the cell centers to appropriate solution variables.
 *
 ********************************************************/
template<class cState, class pState>
void Grid(RadMom1D_UniformMesh<cState, pState> *Soln,
          const double &xMin,
          const double &xMax,
          const int Number_of_Cells) {

  int i;
  double delta_x;

  int TC;

  TC = Number_of_Cells+2*Soln[0].Nghost; // total number of cells

  /* Determine the mesh spacing. */
  delta_x = (xMax - xMin)/double(Number_of_Cells);
  Soln[0].X.setsize(delta_x);

  /* Create the cells. */

  Soln[0].X.x = xMin - (Soln[0].Nghost - HALF)*delta_x;
  // Soln[0].CellHighOrder().AssociateGeometry(Soln[0].X);   // Associate geometry with high-order solution variables


  for ( i = 1 ; i <= TC-1 ; ++i ) {
    // Initialize the coordinate of the centroids
    Soln[i].X.x =  Soln[0].X.x + double(i)*delta_x;
  } /* endfor */
}

#endif /* _RADMOM1D_MESH_INCLUDED  */
