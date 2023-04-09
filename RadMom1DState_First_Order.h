/*!\file RadMom1DState_First_Order.h
  \brief Header file defining 1D RadMom Solution State Classes. */

#ifndef _RADMOM1D_STATE_FIRST_ORDER_INCLUDED
#define _RADMOM1D_STATE_FIRST_ORDER_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

/* Using std namespace functions */
using namespace std;

#ifndef _MATH_MACROS_INCLUDED
#include "./Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _MEDIUM_1D_STATE_INCLUDED
#include "./Medium1DState/Medium1DState.h"
#endif //_MEDIUM_1D_STATE_INCLUDED

/* Define some constants. */
#define STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER 2

// Enviroment flag for cfrte1D root directory path
#define PATHVAR "cfrte1D_Path"

#define TOLER_DENSITY 1e-6
#define TOLER_FLUX 1e-7
#define TOLER_RADIATIVE_DENSITY_M1 1.0e-16

/* Define additional structures. */
struct Eigenstructure_M1 {
    // In order to avoid repeated computations of the eigenvalues whenever we need
    // to compute the eigenvectors, the latter are precomputed and stored in the 
    // current data structure
    double E_dChi2_dE, dChi2_df;
    
    double lambdas[STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER];
    double dFdU[STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER];
    double rc_vec[STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER];
    double lc_vec[STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER];
    
    void Nullify() {
        for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++) {
            lambdas[i] = ZERO;
            for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; j++) {
                dFdU[i][j] = ZERO;
                rc_vec[i][j] = ZERO;
                lc_vec[i][j] = ZERO;
            }
        }
    }
};

/* Define the classes. */
class RadMom1D_cState_First_Order;
class RadMom1D_pState_First_Order;

/*!
 * Class: RadMom1D_pState_First_Order
 *
 * @brief Primitive variable solution state class definition for radiative
 *        transfer.
 *
 * Primitive variable solution state class definition for radiative
 * transfer
 *
 * \verbatim
 * Member functions
 *     e        -- Return radiative energy (density).
 *     f        -- Return normalized radiative flux (vector).
 *     fsca     -- Return absolute normalized flux (scalar).
 *     f2       -- Return square of radiative flux (scalar).
 *     xi       -- Return sqrt(4-3*fsca^2).
 *     chi      -- Return (3+4*fsca^2)/(5+2*xi).
 *     k        -- Return Boltzmann constant.
 *     h        -- Return Planck's constant.
 *     c        -- Return speed of light.
 *     sig      -- Return addition of absorption and scattering coefficients.
 *     a        -- Return 8*pi^5*k^4/(15*h^3*c^3).
 *     U        -- Return conserved solution state.
 *     F        -- Return x-direction solution flux.
 *     Fx       -- Return x-direction solution flux.
 *     Fy       -- Return y-direction solution flux.
 *     dFdU     -- Return x-direction flux Jacobian.
 *     dFxdU    -- Return x-direction flux Jacobian.
 *     dFydU    -- Return y-direction flux Jacobian.
 *     lambda   -- Return x-direction eigenvalue(s).
 *     lambda_x -- Return x-direction eigenvalue(s).
 *     rp       -- Return primitive right eigenvector (x-direction).
 *     rc       -- Return conserved right eigenvector (x-direction).
 *     lp       -- Return primitive left eigenvector (x-direction).
 *     lp_x     -- Return primitive left eigenvector (x-direction).
 *     S        -- Return source term vector.
 *     dUdW     -- Return Jacobian of conserved solution variables wrt
 *                 primitive solution variables.
 *     dWdU     -- Return Jacobian of primitive solution variables wrt
 *                 conserved solution variables.
 *
 * Member operators
 *      W -- a primitive solution state
 *      b -- a scalar (double)
 *
 * W = W;
 * c = W[i];
 * W = W + W;
 * W = W - W;
 * c = W * W; (inner product)
 * W = c * W;
 * W = W * c;
 * W = W / c;
 * W = W ^ W; (a useful product)
 * W = +W;
 * W = -W;
 * W += W;
 * W -= W;
 * W == W;
 * W != W;
 * cout << W; (output function)
 * cin  >> W; (input function)
 * \endverbatim
 */
class RadMom1D_pState_First_Order {
protected:
  
private:
public:
  #ifdef STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER
    double  m_values[STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER];
  #else
    double* m_values;
  #endif
  
  static int              closure_type; //!< Identifier for either M1 or P1 moment closure.
  static int              Absorption_Model;
  static int              Scattering_Func;
  static double                    c; //!< Speed of light.
  static double                    a; //!< Blackbody source constant.
  
  // Planck constants
  static double                   C1; // [W/(m2 ster m-4)]
  static int   NUM_VAR_RADMOM1D_FIRST_ORDER; //!< total number of RadMom1D state variables
  
  //@}
  
  void Setup_Eigenstructure_P1(Eigenstructure_M1 &Eig_P1);
  void Setup_Eigenstructure_P1(Eigenstructure_M1 &Eig_P1) const;
  
  void Setup_Eigenstructure_M1(Eigenstructure_M1 &Eig_M1);
  void Setup_Eigenstructure_M1(Eigenstructure_M1 &Eig_M1) const;
  
  void Setup_Eigenstructure_M1(Eigenstructure_M1 &Eig_M1, const RadMom1D_pState_First_Order &Wl, const RadMom1D_pState_First_Order &Wr);
  void Setup_Eigenstructure_M1(Eigenstructure_M1 &Eig_M1, const RadMom1D_pState_First_Order &Wl, const RadMom1D_pState_First_Order &Wr) const;
  
  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  RadMom1D_pState_First_Order(void) {
      Nullify();
      Allocate();
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
              m_values[i] = ZERO;
      }
  }

  //! Copy constructor.
  RadMom1D_pState_First_Order(const RadMom1D_pState_First_Order &W) {
      Nullify();
      Allocate();
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           m_values[i] = W.m_values[i];  
      }      
  }

  //! Assignment constructor.
  RadMom1D_pState_First_Order(const RadMom1D_cState_First_Order &U);

  //! Assignment constructor.
  RadMom1D_pState_First_Order(const double &den,
                              const double &fx) {
      Nullify();
      Allocate();

      m_values[0] = den;
      m_values[1] = fx; 
  }
  
  //! Value Constructor
  explicit RadMom1D_pState_First_Order(const double &Val);

  /* Destructor. */
  ~RadMom1D_pState_First_Order(void) { Deallocate(); }
  // Use automatically generated destructor.
  //@}

  //@{ @name Useful operators.
  //! Return the number of variables.
  static int StaticNumVar(void) { return STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; }
  static int NumVar(void) { return NUM_VAR_RADMOM1D_FIRST_ORDER; }
  
  //! Copy operator.
  void Copy(const RadMom1D_pState_First_Order &W) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           m_values[i] = W.m_values[i];  
      }
  }
  
  void Copy_to_W(RadMom1D_pState_First_Order &W) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           W.m_values[i] = m_values[i];
      }
  }
  
  void Copy_to_W(RadMom1D_pState_First_Order &W) const {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           W.m_values[i] = m_values[i];
      }
  }

  //! Set array pointers to null
  void Nullify();
  
  //! Vacuum operator.
  void Vacuum(void) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           m_values[i] = ZERO;  
      }      
  }
  //! Ones operator.
  void Ones(void) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           m_values[i] = ONE;
      }
  }
  
  double I0(void) {return m_values[0];}
  double I0(void) const {return m_values[0];}
  
  double N1x(void) {return m_values[1];}
  double N1x(void) const {return m_values[1];}
  
  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const {
      if (closure_type == MOMENT_CLOSURE_M1) {
          if (I0() < ZERO || fsca() > ONE) return 1;
      }
      return 0;
  }
  
  //@}

  //@{ @name State functions.
  //! Absolute flux (scalar).
  double fsca(void);
  double fsca(void) const;

  //! Radiative flux squared (scalar).
  double f2(void);
  double f2(void) const;

 //! xi.
  double xi(void);
  double xi(void) const;

  //! chi.
  double chi(void);
  double chi(void) const;
  
  //@{ Routines for computations of the interpolative-based approximations of the 
  // Eddingtong factor for the non-gray M1 closure obeying Bose-Einstein statistics
  void Evaluate_Chi2_derivatives_Gray(double &E_dChi2_dE, double &dChi2_df);
  void Evaluate_Chi2_derivatives_Gray(double &E_dChi2_dE, double &dChi2_df) const;
  
  void Evaluate_Chi2_derivatives(double &E_dChi2_dE, double &dChi2_df);
  void Evaluate_Chi2_derivatives(double &E_dChi2_dE, double &dChi2_df) const;
  //@} 
  
  double Evaluate_Lambdas_x(const int &index_U);
  double Evaluate_Lambdas_x(const int &index_U) const;
  void Evaluate_Lambdas_x(RadMom1D_pState_First_Order &lambda);
  void Evaluate_Lambdas_x(RadMom1D_pState_First_Order &lambda) const;
  
  void Evaluate_Lambdas_x(Eigenstructure_M1 &Eig_M1);
  void Evaluate_Lambdas_x(Eigenstructure_M1 &Eig_M1) const;
  
  //@{ @name Conserved solution state.
  double U(const int &index_U);
  double U(const int &index_U) const;
  RadMom1D_cState_First_Order U(void);
  RadMom1D_cState_First_Order U(void) const;
  RadMom1D_cState_First_Order U(const RadMom1D_pState_First_Order &W);
  friend RadMom1D_cState_First_Order U(const RadMom1D_pState_First_Order &W);
  //@}

  //@{ @name Solution flux and Jacobian (x-direction).
  RadMom1D_cState_First_Order Fx(void);
  RadMom1D_cState_First_Order Fx(void) const;
  RadMom1D_cState_First_Order Fx(const RadMom1D_pState_First_Order &W);
  friend RadMom1D_cState_First_Order Fx(const RadMom1D_pState_First_Order &W);
  void dFxdU(Eigenstructure_M1 &Eig_M1);
  void dFxdU(Eigenstructure_M1 &Eig_M1) const;
  //@}
  
  //@{ @name Eigenvalue(s) (x-direction).
  RadMom1D_pState_First_Order lambda_x(void);
  RadMom1D_pState_First_Order lambda_x(void) const;
  RadMom1D_pState_First_Order lambda_x(const RadMom1D_pState_First_Order &W);
  friend RadMom1D_pState_First_Order lambda_x(const RadMom1D_pState_First_Order &W);
  double lambda_x(int index);
  double lambda_x(int index) const;
  friend double lambda_x(const RadMom1D_pState_First_Order &W, int index);
  RadMom1D_pState_First_Order lambda_x(const double &V) const;
  //@}

  void Compute_Correction_Factor_Roe(double *Correction_Factor,
                                     const RadMom1D_pState_First_Order &Wl,
                                     const RadMom1D_pState_First_Order &Wr);
  
  void Compute_Correction_Factor_Roe(double *Correction_Factor,
                                     const RadMom1D_pState_First_Order &Wl,
                                     const RadMom1D_pState_First_Order &Wr) const;
                                     
  void Compute_Correction_Factor_Approximate_Jacobian_Roe(double A_Roe[][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER],
                                                          const RadMom1D_pState_First_Order &Wl,
                                                          const RadMom1D_pState_First_Order &Wr);
                                     
  void Compute_Correction_Factor_Approximate_Jacobian_Roe(double A_Roe[][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER],
                                                          const RadMom1D_pState_First_Order &Wl,
                                                          const RadMom1D_pState_First_Order &Wr) const;

  //@{ @name Conserved right eigenvector (x-direction).
  void rc(Eigenstructure_M1 &Eig_M1);
  void rc(Eigenstructure_M1 &Eig_M1) const;
  friend void rc(const RadMom1D_pState_First_Order &W, Eigenstructure_M1 &Eig_M1);

  //@{ @name Conserved left eigenvector (x-direction).
  void lc(Eigenstructure_M1 &Eig_M1);
  void lc(Eigenstructure_M1 &Eig_M1) const;
  friend void lc(const RadMom1D_pState_First_Order &W, Eigenstructure_M1 &Eig_M1);
  //@}

  //@{ @name First order linear solution reconstruction
  void Reconstruct( const RadMom1D_pState_First_Order &Wc,
                    const RadMom1D_pState_First_Order &phi,
                    const RadMom1D_pState_First_Order &dWdx,
                    const double &dX);
  
    //! @name Allocators and deallocators
  //@{

  //! memory allocation / deallocation for the I array
  void Allocate();
  void Deallocate();
  
  //! Setup state static variables
  //! memory allocation / deallocation for the static arrays
  static void SetupStatic( const int &ScatteringFunc,
                           const int &i_Moment_Closure,
                           const int &i_AbsorptionModel);
  static void DeallocateStatic();
  
  
  void RoeAverage(const RadMom1D_pState_First_Order &Wl,
                  const RadMom1D_pState_First_Order &Wr);
  
  void AverageStates(const RadMom1D_pState_First_Order &Wl,
                     const RadMom1D_pState_First_Order &Wr);
  
  void Rotate(const double &norm_dir_x);
  
  void Characteristic(RadMom1D_pState_First_Order W_wall,
                      const double &wall_temperature, 
                      const double &wall_emissivity, 
                      const double &norm_dir);
  void Gray_Wall(RadMom1D_pState_First_Order W_wall,
                 const double &wall_temperature, 
                 const double &wall_emissivity, 
                 const double &norm_dir);

  void Partial_Normalized_Moments_1D(double &E_plus, double &Fx_plus);

  void PartialMoments_n(RadMom1D_pState_First_Order W_inner,
                        const double &wall_temperature, 
                        const double &wall_emissivity, 
                        const double &norm_dir);
  
  //! Index operator.
  double &operator[](int index)       { 
    assert( index >= 1 && index <= NumVar() );
    return m_values[index-1]; 
  }
  const double &operator[](int index) const { 
    assert( index >= 1 && index <= NumVar() );
    return m_values[index-1]; 
  }
  
  //@{ @name Binary arithmetic operators.
  RadMom1D_pState_First_Order operator +(const RadMom1D_pState_First_Order &W) const;
  RadMom1D_pState_First_Order operator -(const RadMom1D_pState_First_Order &W) const;
  RadMom1D_pState_First_Order operator *(const double &b) const;
  double operator *(const RadMom1D_pState_First_Order &W) const;
  friend RadMom1D_pState_First_Order operator *(const double &b, const RadMom1D_pState_First_Order &W);
  RadMom1D_pState_First_Order operator /(const double &b) const;
  RadMom1D_pState_First_Order operator /(const RadMom1D_pState_First_Order &W) const;
  RadMom1D_pState_First_Order operator ^(const RadMom1D_pState_First_Order &W) const;
  
  /* Assignment operator. */
  RadMom1D_pState_First_Order &operator = (const RadMom1D_pState_First_Order &W);
  
  //@{ @name Shortcut arithmetic operators.
  RadMom1D_pState_First_Order &operator +=(const RadMom1D_pState_First_Order &W);
  RadMom1D_pState_First_Order &operator -=(const RadMom1D_pState_First_Order &W);
  RadMom1D_pState_First_Order &operator *=(const double &b);
  RadMom1D_pState_First_Order &operator *=(const RadMom1D_pState_First_Order &W);
  RadMom1D_pState_First_Order &operator /=(const double &b);
  RadMom1D_pState_First_Order &operator /=(const RadMom1D_pState_First_Order &W);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const RadMom1D_pState_First_Order &W);
  friend istream &operator >> (istream &in_file,  RadMom1D_pState_First_Order &W);
  //@}

  // Routines for computing approximate Roe matrix
  RadMom1D_pState_First_Order U_to_W_Roe(const RadMom1D_cState_First_Order &U);
  RadMom1D_pState_First_Order U_to_W_Roe(const RadMom1D_cState_First_Order &U) const;
  
  RadMom1D_cState_First_Order W_Roe_to_U(const RadMom1D_pState_First_Order &W_Roe);
  RadMom1D_cState_First_Order W_Roe_to_U(const RadMom1D_pState_First_Order &W_Roe) const;
  
  void Approximate_Roe_Matrix(Eigenstructure_M1 &Eig_M1,
                              const RadMom1D_pState_First_Order &Wl,
                              const RadMom1D_pState_First_Order &Wr);
  
  void Approximate_Roe_Matrix(Eigenstructure_M1 &Eig_M1,
                              const RadMom1D_pState_First_Order &Wl,
                              const RadMom1D_pState_First_Order &Wr) const;

};

 // end of RadMom1D_pState_First_Order class

/*!
 * Class: RadMom1D_cState_First_Order
 *
 * @brief Conserved variable solution state class definition for radiative
 *        transfer.
 *
 * Conserved variable solution state class definition for radiative
 * transfer.
 *
 * \verbatim
 * Member functions
 *     E        -- Return radiative energy (density).
 *     F        -- Return radiative flux (vector).
 *     absorp   -- Return absorption coefficient.
 *     emiss    -- Return emission coefficient.
 *     scatt    -- Return scattering coefficient.
 *     fsca     -- Return absolute normalized flux (scalar).
 *     f2       -- Return square of radiative flux (scalar).
 *     xi       -- Return sqrt(4-3*fsca^2).
 *     chi      -- Return (3+4*fsca^2)/(5+2*xi).
 *     k        -- Return Boltzmann constant.
 *     h        -- Return Planck's constant.
 *     c        -- Return speed of light.
 *     sig      -- Return addition of absorption and scattering coefficients.
 *     a        -- Return 8*pi^5*k^4/(15*h^3*c^3).
 *     W        -- Return primitive solution state.
 *     Flux     -- Return x-direction solution flux.
 *     Fx       -- Return x-direction solution flux.
 *     Fy       -- Return y-direction solution flux.
 *     dFdU     -- Return x-direction flux Jacobian.
 *     dFxdU    -- Return x-direction flux Jacobian.
 *     dFydU    -- Return y-direction flux Jacobian.
 *     S        -- Return axisymmetric source term vector.
 *     dUdW     -- Return Jacobian of conserved solution variables wrt
 *                 primitive solution variables.
 *     dWdU     -- Return Jacobian of primitive solution variables wrt
 *                 conserved solution variables.
 * Member operators
 *      U -- a primitive solution state
 *      b -- a scalar (double)
 *
 * U = U;
 * c = U[i];
 * U = U + U;
 * U = U - U;
 * c = U * U; (inner product)
 * U = c * U;
 * U = U * c;
 * U = U / c;
 * U = U ^ U; (a useful product)
 * U = +U;
 * U = -U;
 * U += U;
 * U -= U;
 * U == U;
 * U != U;
 * cout << U; (output function)
 * cin  >> U; (input function)
 * \endverbatim
 */
class RadMom1D_cState_First_Order {
protected:
private:
public:
  
  #ifdef STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER
  double  m_values[STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER];
  #else 
  double* m_values;
  #endif
  
  static int              closure_type; //!< Identifier for either M1 or P1 moment closure.
  static int          Absorption_Model;
  static int          Scattering_Func;
  static double                      c; //!< Speed of light.
  static double                      a; //!< Blackbody source constant.
  
  // Planck constants
  static double                     C1; // [W/(m2 ster m-4)]static int   NUM_VAR_RADMOM1D_FIRST_ORDER; //!< total number of RadMom1D state variables

  static int   NUM_VAR_RADMOM1D_FIRST_ORDER; //!< total number of RadMom1D state variables

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  RadMom1D_cState_First_Order(void) {
      Nullify();
      Allocate();
       
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           m_values[i] = ZERO;
      }
  }

  //! Copy constructor.
  RadMom1D_cState_First_Order(const RadMom1D_cState_First_Order &U) {
      Nullify();
      Allocate();
      
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           m_values[i] = U.m_values[i];  
      }
  }

  //! Copy constructor.
  RadMom1D_cState_First_Order(const RadMom1D_pState_First_Order &W);

  //! Value Constructor
  explicit RadMom1D_cState_First_Order(const double &Val);

  //! Assignment constructor.
  RadMom1D_cState_First_Order(const double &den,
		  const double &flux) {
      Nullify();
      Allocate();
      
      m_values[0] = den;
      m_values[1] = flux; 
  }
  
  /* Destructor. */
  ~RadMom1D_cState_First_Order(void) { Deallocate(); }
  // Use automatically generated destructor.
  //@}

  //@{ @name Useful operators.
  //! Return the number of variables.
  static int StaticNumVar(void) { return STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; }
  static int NumVar(void) { return NUM_VAR_RADMOM1D_FIRST_ORDER; }
  
  //! Copy operator.
  void Copy(const RadMom1D_cState_First_Order &U) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
          m_values[i] = U.m_values[i];  
      }
  }
  
  void Copy_to_U(RadMom1D_cState_First_Order &U) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           U.m_values[i] = m_values[i];
      }
  }
  
  void Copy_to_U(RadMom1D_cState_First_Order &U) const {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           U.m_values[i] = m_values[i];
      }
  }
  
  //! Set array pointers to null
  void Nullify();
  
  //! Vacuum operator.
  void Vacuum(void) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           m_values[i] = ZERO;  
      }
  }
  //! Ones operator.
  void Ones(void) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
           m_values[i] = ONE;
      }
  }
  
  double I0(void) {return m_values[0];}
  double I0(void) const {return m_values[0];}
  
  double I1x(void) {return m_values[1];}
  double I1x(void) const {return m_values[1];}
  
  //! Absolute radiative flux (scalar).
  double fsca(void);
  double fsca(void) const;

  //! Square of radiative flux (scalar).
  double f2(void);
  double f2(void) const;
  
  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const {
      return W().Unphysical_Properties();
  }    

  //@{ @name State functions.
  //! xi.
  double xi(void);
  double xi(void) const;

  //! chi.
  double chi(void);
  double chi(void) const;
  
    //@{ @name Primitive solution state.
  RadMom1D_pState_First_Order W(void);
  RadMom1D_pState_First_Order W(void) const;
  RadMom1D_pState_First_Order W(const RadMom1D_cState_First_Order &U);
  friend RadMom1D_pState_First_Order W(const RadMom1D_cState_First_Order &U);
  //@}
  
  //@{ @name Solution flux and Jacobian (x-direction).
  RadMom1D_cState_First_Order Fx(void);
  RadMom1D_cState_First_Order Fx(void) const;
  RadMom1D_cState_First_Order Fx(const RadMom1D_cState_First_Order &U);
  friend RadMom1D_cState_First_Order Fx(const RadMom1D_cState_First_Order &U);
  //@}

  //@{ @name Eigenvalue(s) (x-direction).
  double lambda(const int &index);
  double lambda(const int &index) const;
  friend double lambda(const RadMom1D_cState_First_Order &U, const int &index);
  
  //@{ @name Include all source vectors and Jacobians.
  RadMom1D_cState_First_Order S(const Medium1D_State &M);
  double Sr(const Medium1D_State &M);
  double Sr(const Medium1D_State &M) const; // add const declaration and definition for outputting in RadMom1D_Mesh.h

  //! memory allocation / deallocation for the I array
  void Allocate();
  void Deallocate();
  
  //! Setup state static variables
  //! memory allocation / deallocation for the static arrays
  static void SetupStatic( const int &ScatteringFunc,
                           const int &i_Moment_Closure,
                           const int &i_AbsorptionModel);
  static void DeallocateStatic();

  void Rotate(const double &norm_dir_x);
  
  void Set_ICs(const double &Medium_Temperature);
  void Set_ICs_Beam(const double &Medium_Temperature, const int &Direction);
  void Set_ICs_Intensity(const double &Ib_wall);
  void Set_BCs(const double *Intensity, const double norm_dir);

  //! Index operator.
  double &operator[](int index) {
    assert( index >= 1 && index <= NumVar() );
    return m_values[index-1]; 
  }
  const double &operator[](int index) const {
    assert( index >= 1 && index <= NumVar() );
    return m_values[index-1]; 
  }
  
  //@{ @name Binary arithmetic operators.
  RadMom1D_cState_First_Order operator +(const RadMom1D_cState_First_Order &W) const;
  RadMom1D_cState_First_Order operator -(const RadMom1D_cState_First_Order &W) const;
  RadMom1D_cState_First_Order operator *(const double &b) const;
  double operator *(const RadMom1D_cState_First_Order &W) const;
  friend RadMom1D_cState_First_Order operator *(const double &b, const RadMom1D_cState_First_Order &W);
  RadMom1D_cState_First_Order operator /(const double &b) const;
  RadMom1D_cState_First_Order operator /(const RadMom1D_cState_First_Order &W) const;
  RadMom1D_cState_First_Order operator ^(const RadMom1D_cState_First_Order &W) const;
  
  /* Assignment operator. */
  RadMom1D_cState_First_Order &operator = (const RadMom1D_cState_First_Order &W);
  
  //@{ @name Shortcut arithmetic operators.
  RadMom1D_cState_First_Order &operator +=(const RadMom1D_cState_First_Order &W);
  RadMom1D_cState_First_Order &operator -=(const RadMom1D_cState_First_Order &W);
  RadMom1D_cState_First_Order &operator *=(const double &b);
  RadMom1D_cState_First_Order &operator *=(const RadMom1D_cState_First_Order &W);
  RadMom1D_cState_First_Order &operator /=(const double &b);
  RadMom1D_cState_First_Order &operator /=(const RadMom1D_cState_First_Order &W);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const RadMom1D_cState_First_Order &U);
  friend istream &operator >> (istream &in_file,  RadMom1D_cState_First_Order &U);
  //@}
};

// End of RadMom1D_cState_First_Order class

inline void RadMom1D_pState_First_Order :: Evaluate_Lambdas_x(RadMom1D_pState_First_Order &lambda) {
    double wave;
    double E_dChi2_dE, dChi2_df;
    
    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            lambda[1] = -sqrt(THREE)/THREE;
            lambda[2] = sqrt(THREE)/THREE;
            break;
        case MOMENT_CLOSURE_M1 :
            // if (f2() < TOLER_FLUX){
            //     lambda[1] = -sqrt(THREE)/THREE;
            //     lambda[2] = sqrt(THREE)/THREE;
            // } else {
                // Compute derivatives of the second-order closing fluxes with respect to the
                // lower-order angular moments
                Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_df);
                wave = sqrt(sqr(HALF*dChi2_df) + chi() - N1x()*dChi2_df);
                lambda[1] = HALF*dChi2_df - wave;
                lambda[2] = HALF*dChi2_df + wave;
            // }
            break;
    };
}

inline void RadMom1D_pState_First_Order :: Evaluate_Lambdas_x(RadMom1D_pState_First_Order &lambda) const {
    double wave;
    double E_dChi2_dE, dChi2_df;

    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            lambda[1] = -sqrt(THREE)/THREE;
            lambda[2] = sqrt(THREE)/THREE;
            break;
        case MOMENT_CLOSURE_M1 :
            // if (f2() < TOLER_FLUX){
            //     lambda[1] = -sqrt(THREE)/THREE;
            //     lambda[2] = sqrt(THREE)/THREE;
            // } else {
                // Compute derivatives of the second-order closing fluxes with respect to the
                // lower-order angular moments
                Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_df);
                wave = sqrt(sqr(HALF*dChi2_df) + chi() - N1x()*dChi2_df);
                lambda[1] = HALF*dChi2_df - wave;
                lambda[2] = HALF*dChi2_df + wave;
            // }
            break;
    };
}

inline void RadMom1D_pState_First_Order :: Evaluate_Lambdas_x(Eigenstructure_M1 &Eig_M1) {
    double wave;
    double E_dChi2_dE, dChi2_df;
    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            Eig_M1.lambdas[0] = -sqrt(THREE)/THREE;
            Eig_M1.lambdas[1] = sqrt(THREE)/THREE;
            break;
        case MOMENT_CLOSURE_M1 :
            // if (f2() < TOLER_FLUX){
            //     Eig_M1.lambdas[0] = -sqrt(THREE)/THREE;
            //     Eig_M1.lambdas[1] = sqrt(THREE)/THREE;
            // } else {
                // Compute derivatives of the second-order closing fluxes with respect to the
                // lower-order angular moments
                Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_df);
                wave = sqrt(sqr(HALF*dChi2_df) + chi() - N1x()*dChi2_df);
                Eig_M1.lambdas[0] = HALF*dChi2_df - wave;
                Eig_M1.lambdas[1] = HALF*dChi2_df + wave;
            // }
            break;
        default:
            cout << "Closure type not specified" << endl;
            exit(0);
            break;
    }
}

inline void RadMom1D_pState_First_Order :: Evaluate_Lambdas_x(Eigenstructure_M1 &Eig_M1) const {
    double wave;
    double E_dChi2_dE, dChi2_df;
    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            Eig_M1.lambdas[0] = -sqrt(THREE)/THREE;
            Eig_M1.lambdas[1] = sqrt(THREE)/THREE;
            break;
        case MOMENT_CLOSURE_M1 :
            // if (f2() < TOLER_FLUX){
            //     Eig_M1.lambdas[0] = -sqrt(THREE)/THREE;
            //     Eig_M1.lambdas[1] = sqrt(THREE)/THREE;
            // } else {
                // Compute derivatives of the second-order closing fluxes with respect to the
                // lower-order angular moments
                Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_df);
                wave = sqrt(sqr(HALF*dChi2_df) + chi() - N1x()*dChi2_df);
                Eig_M1.lambdas[0] = HALF*dChi2_df - wave;
                Eig_M1.lambdas[1] = HALF*dChi2_df + wave;
            // }
            break;
        default:
            cout << "Closure type not specified" << endl;
            exit(0);
            break;
    }
}

//*********************************************************************************
// Routines for the interpolative-based approximation of the non-gray M1 closure.
// ********************************************************************************
inline void RadMom1D_pState_First_Order::Evaluate_Chi2_derivatives_Gray(double &E_dChi2_dE, double &dChi2_df) {
    double norm_f = fsca();
    dChi2_df = TWO * norm_f / xi();
    E_dChi2_dE = -norm_f*dChi2_df;
}

inline void RadMom1D_pState_First_Order::Evaluate_Chi2_derivatives_Gray(double &E_dChi2_dE, double &dChi2_df) const {
    double norm_f = fsca();
    dChi2_df = TWO * norm_f / xi();
    E_dChi2_dE = -norm_f*dChi2_df;
}

inline void RadMom1D_pState_First_Order::Evaluate_Chi2_derivatives(double &E_dChi2_dE, double &dChi2_df) {
    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            E_dChi2_dE = ZERO;
            dChi2_df = ZERO;
            break;
        case MOMENT_CLOSURE_M1 :
            Evaluate_Chi2_derivatives_Gray(E_dChi2_dE, dChi2_df);
            break; 
    }
}

inline void RadMom1D_pState_First_Order::Evaluate_Chi2_derivatives(double &E_dChi2_dE, double &dChi2_df) const {
    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            E_dChi2_dE = ZERO;
            dChi2_df = ZERO;
            break;
        case MOMENT_CLOSURE_M1 :
            Evaluate_Chi2_derivatives_Gray(E_dChi2_dE, dChi2_df);
            break; 
    }
}

/****************************************************************//**
 * Set array pointers to null.
 ********************************************************************/
inline void RadMom1D_pState_First_Order :: Nullify() {
    #ifndef STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER
    m_values = NULL;
    #endif
};

/****************************************************************//**
 * Array allocator and deallocator for intensity array.
 ********************************************************************/
inline void RadMom1D_pState_First_Order :: Allocate() {
    #ifndef STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER
    //   deallocate first
    Deallocate();

    //   create the jagged array
    if (NUM_VAR_RADMOM1D_FIRST_ORDER > 0)  {
        m_values = new double[NUM_VAR_RADMOM1D_FIRST_ORDER];
    }
    #endif
}


inline void RadMom1D_pState_First_Order :: Deallocate() {
    #ifndef STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER
    if ( m_values != NULL ) {
        delete[] m_values; m_values = NULL;
    }
    #endif
}

inline void RadMom1D_pState_First_Order :: SetupStatic( const int &ScatteringFunc,
                                                        const int &i_Moment_Closure,
                                                        const int &i_AbsorptionModel)
{
    // deallocate static vars first, just to be safe
    DeallocateStatic();
    
    NUM_VAR_RADMOM1D_FIRST_ORDER = STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER;
    
    closure_type = i_Moment_Closure;
    Absorption_Model = i_AbsorptionModel;
    Scattering_Func = ScatteringFunc;

    if (closure_type != MOMENT_CLOSURE_P1 &&
        closure_type != MOMENT_CLOSURE_M1) {
        cout << "Closure type not specified properly !!!!" << endl;
        cout << "Should specify MOMENT_CLOSURE_M1 or MOMENT_CLOSURE_P1" << endl;
        exit(0);
    }
}

inline void RadMom1D_pState_First_Order :: DeallocateStatic()
{

}

/********************************************
 * RadMom1D_pState_First_Order Value Constructor.        *
 *******************************************/
inline RadMom1D_pState_First_Order::RadMom1D_pState_First_Order(const double &Val){
   for ( int i=0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
     m_values[i] = Val;
  }
}


/********************************************************
 * RadMom1D_pState_First_Order::fsca -- Absolute radiative flux.    *
 ********************************************************/
inline double RadMom1D_pState_First_Order::fsca() {
  return ( sqrt(f2()) );
}

inline double RadMom1D_pState_First_Order::fsca() const {
  return ( sqrt(f2()) );
}

/********************************************************
 * RadMom1D_pState_First_Order::f2 -- Square of radiative flux.     *
 ********************************************************/
inline double RadMom1D_pState_First_Order::f2() {
    double f2_val;
    f2_val = sqr(N1x());
    return f2_val;
}

inline double RadMom1D_pState_First_Order::f2() const {
    double f2_val;
    f2_val = sqr(N1x());
    return f2_val;
}

/********************************************************
 * RadMom1D_pState_First_Order::xi -- Useful variable.              *
 ********************************************************/
inline double RadMom1D_pState_First_Order::xi() {
  return (sqrt(FOUR-THREE*f2()));
}

inline double RadMom1D_pState_First_Order::xi() const {
  return (sqrt(FOUR-THREE*f2()));
}

/********************************************************
 * RadMom1D_pState_First_Order::chi -- Useful ratio.                *
 ********************************************************/
inline double RadMom1D_pState_First_Order::chi() {
    double Chi2;
    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            Chi2 = ONE/THREE;
            break;
        case MOMENT_CLOSURE_M1 :
            Chi2 = (THREE+FOUR*f2())/(FIVE+TWO*xi());
            break;
        default:
            cout << "Closure type not specified" << endl;
            exit(0);
            break;
    };
    return (Chi2);
}

inline double RadMom1D_pState_First_Order::chi() const {
    double Chi2;
    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            Chi2 = ONE/THREE;
            break;
        case MOMENT_CLOSURE_M1 :
            Chi2 = (THREE+FOUR*f2())/(FIVE+TWO*xi());
            break;
        default:
            cout << "Closure type not specified" << endl;
            exit(0);
            break;
    };
    return (Chi2);
}

/****************************************************************//**
 * RadMom1D_pState_First_Order -- Binary arithmetic operators.
 ********************************************************************/
//! addition
inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order :: operator +(const RadMom1D_pState_First_Order &W) const {
  RadMom1D_pState_First_Order W_temp(*this);
  W_temp += W;
  return W_temp;
}

//! subtraction
inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order :: operator -(const RadMom1D_pState_First_Order &W) const {
  RadMom1D_pState_First_Order W_temp(*this);
  W_temp -= W;
  return W_temp;
}

//! scalar multiplication
inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order :: operator *(const double &b) const {
  RadMom1D_pState_First_Order W_temp(*this);
  W_temp *= b;
  return W_temp;
}

inline RadMom1D_pState_First_Order operator *(const double &b, const RadMom1D_pState_First_Order &W) {
  RadMom1D_pState_First_Order W_temp(W);
  W_temp *= b;
  return W_temp;
}

//! scalar division
inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order :: operator /(const double &b) const {
  RadMom1D_pState_First_Order W_temp(*this);
  W_temp /= b;
  return W_temp;
}

//! solution state division operator
inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order :: operator /(const RadMom1D_pState_First_Order &W) const {
  RadMom1D_pState_First_Order W_temp(*this);
  W_temp /= W;
  return W_temp;
}

//! inner product
inline double RadMom1D_pState_First_Order :: operator *(const RadMom1D_pState_First_Order &W) const {
  double sum=0.0;

  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    sum += m_values[i]*W.m_values[i];
  }
  return sum;
}

//! solution state product operator
inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order :: operator ^(const RadMom1D_pState_First_Order &W) const {
  RadMom1D_pState_First_Order W_temp(*this);
  W_temp *= W;
  return W_temp;
}

/****************************************************************//**
 * RadMom1D_pState_First_Order -- Assignment operator.
 ********************************************************************/
inline RadMom1D_pState_First_Order& RadMom1D_pState_First_Order :: operator =(const RadMom1D_pState_First_Order &W) {
  if( this != &W) Copy(W);
  return (*this);
}

/****************************************************************//**
 * RadMom1D_pState_First_Order -- Shortcut operators.
 ********************************************************************/
inline RadMom1D_pState_First_Order& RadMom1D_pState_First_Order :: operator +=(const RadMom1D_pState_First_Order &W) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] += W.m_values[i];
  }
  return (*this);
}

inline RadMom1D_pState_First_Order& RadMom1D_pState_First_Order :: operator -=(const RadMom1D_pState_First_Order &W) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] -= W.m_values[i];
  }
    return (*this);
}

inline RadMom1D_pState_First_Order& RadMom1D_pState_First_Order::operator *=(const double &b) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] *= b;
  }
    return (*this);
}

inline RadMom1D_pState_First_Order& RadMom1D_pState_First_Order::operator *=(const RadMom1D_pState_First_Order &W) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] *= W.m_values[i];
  }
    return (*this);
}

inline RadMom1D_pState_First_Order& RadMom1D_pState_First_Order::operator /=(const double &b) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] /= b;
  }
  return (*this);
}

inline RadMom1D_pState_First_Order& RadMom1D_pState_First_Order::operator /=(const RadMom1D_pState_First_Order &W) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] /= W.m_values[i];
  }
  return (*this);
}

/********************************************************
 * RadMom1D_pState_First_Order -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const RadMom1D_pState_First_Order &W) {
    out_file.setf(ios::scientific);

    out_file << " " << W[1] << " " << W[2] << "\n";

    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream &operator >> (istream &in_file, RadMom1D_pState_First_Order &W) {
    in_file.setf(ios::skipws);

    in_file >> W[1] >> W[2];

    in_file.unsetf(ios::skipws);
    return (in_file);
}

/****************************************************************//**
 * Set array pointers to null.
 ********************************************************************/
inline void RadMom1D_cState_First_Order :: Nullify() 
{ 
#ifndef STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER
  m_values = NULL;
#endif
};

/****************************************************************//**
 * Array allocator and deallocator for intensity array.
 ********************************************************************/
inline void RadMom1D_cState_First_Order :: Allocate()
{
#ifndef STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER

  // deallocate first
  Deallocate();

  // create the jagged array
  if (NUM_VAR_RADMOM1D_FIRST_ORDER > 0)  { 
      m_values = new double[NUM_VAR_RADMOM1D_FIRST_ORDER]; 
  }

#endif
}


inline void RadMom1D_cState_First_Order :: Deallocate()
{

#ifndef STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER
  if ( m_values != NULL ) { 
      delete[] m_values; m_values = NULL;
  }
#endif
}

inline void RadMom1D_cState_First_Order :: SetupStatic( const int &ScatteringFunc,
                                                        const int &i_Moment_Closure,
                                                        const int &i_AbsorptionModel)
{
    // deallocate static vars first, just to be safe
    DeallocateStatic();
    
    NUM_VAR_RADMOM1D_FIRST_ORDER = STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER;

    closure_type = i_Moment_Closure;
    Absorption_Model = i_AbsorptionModel;
    Scattering_Func = ScatteringFunc;

    if (closure_type != MOMENT_CLOSURE_P1 &&
        closure_type != MOMENT_CLOSURE_M1) {
        cout << "Closure type not specified properly !!!!" << endl;
        cout << "Should specify MOMENT_CLOSURE_M1 or MOMENT_CLOSURE_P1" << endl;
        exit(0);
    }
}

inline void RadMom1D_cState_First_Order :: DeallocateStatic()
{

}

/********************************************
 * RadMom1D_cState_First_Order Value Constructor.        *
 *******************************************/
inline RadMom1D_cState_First_Order::RadMom1D_cState_First_Order(const double &Val){
    for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ){
      m_values[i] = Val;
    }
}

/********************************************************
 * RadMom1D_cState_First_Order::fsca -- Absolute radiative flux.    *
 ********************************************************/
inline double RadMom1D_cState_First_Order::fsca() {
  return ( sqrt(f2()) );
}

inline double RadMom1D_cState_First_Order::fsca() const {
  return ( sqrt(f2()) );
}

/********************************************************
 * RadMom1D_cState_First_Order::f2 -- Square of radiative flux.     *
 ********************************************************/
inline double RadMom1D_cState_First_Order::f2() {
  return W().f2();
}

inline double RadMom1D_cState_First_Order::f2() const {
  return W().f2();
}

/********************************************************
 * RadMom1D_cState_First_Order::xi -- Useful variable.              *
 ********************************************************/
inline double RadMom1D_cState_First_Order::xi() {
  return (sqrt(FOUR-THREE*f2()));
}

inline double RadMom1D_cState_First_Order::xi() const {
  return (sqrt(FOUR-THREE*f2()));
}

/********************************************************
 * RadMom1D_cState_First_Order::chi -- Useful ratio.                *
 ********************************************************/
inline double RadMom1D_cState_First_Order::chi() {
    return W().chi();
}

inline double RadMom1D_cState_First_Order::chi() const {
    return W().chi();
}

/****************************************************************//**
 * RadMom1D_cState_First_Order -- Binary arithmetic operators.
 ********************************************************************/
//! addition
inline RadMom1D_cState_First_Order RadMom1D_cState_First_Order :: operator +(const RadMom1D_cState_First_Order &U) const {
  RadMom1D_cState_First_Order U_temp(*this);
  U_temp += U;
  return U_temp;
}

//! subtraction
inline RadMom1D_cState_First_Order RadMom1D_cState_First_Order :: operator -(const RadMom1D_cState_First_Order &U) const {
  RadMom1D_cState_First_Order U_temp(*this);
  U_temp -= U;
  return U_temp;
}

//! scalar multiplication
inline RadMom1D_cState_First_Order RadMom1D_cState_First_Order :: operator *(const double &b) const {
  RadMom1D_cState_First_Order U_temp(*this);
  U_temp *= b;
  return U_temp;
}

inline RadMom1D_cState_First_Order operator *(const double &b, const RadMom1D_cState_First_Order &U) {
  RadMom1D_cState_First_Order U_temp(U);
  U_temp *= b;
  return U_temp;
}

//! scalar division
inline RadMom1D_cState_First_Order RadMom1D_cState_First_Order :: operator /(const double &b) const {
  RadMom1D_cState_First_Order U_temp(*this);
  U_temp /= b;
  return U_temp;
}

//! solution state division operator
inline RadMom1D_cState_First_Order RadMom1D_cState_First_Order :: operator /(const RadMom1D_cState_First_Order &U) const {
  RadMom1D_cState_First_Order U_temp(*this);
  U_temp /= U;
  return U_temp;
}

//! inner product
inline double RadMom1D_cState_First_Order :: operator *(const RadMom1D_cState_First_Order &U) const {
    double sum=0.0;

    for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
      sum += m_values[i]*U.m_values[i];
    }
    return sum;
}

//! solution state product operator
inline RadMom1D_cState_First_Order RadMom1D_cState_First_Order :: operator ^(const RadMom1D_cState_First_Order &U) const {
  RadMom1D_cState_First_Order U_temp(*this);
  U_temp *= U;
  return U_temp;
}

/****************************************************************//**
 * RadMom1D_cState_First_Order -- Assignment operator.
 ********************************************************************/
inline RadMom1D_cState_First_Order& RadMom1D_cState_First_Order :: operator =(const RadMom1D_cState_First_Order &U) {
  if( this != &U) Copy(U);
  return (*this);
}

/****************************************************************//**
 * RadMom1D_cState_First_Order -- Shortcut operators.
 ********************************************************************/
inline RadMom1D_cState_First_Order& RadMom1D_cState_First_Order :: operator +=(const RadMom1D_cState_First_Order &U) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] += U.m_values[i];
  }
  return (*this);
}

inline RadMom1D_cState_First_Order& RadMom1D_cState_First_Order :: operator -=(const RadMom1D_cState_First_Order &U) {
    for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
      m_values[i] -= U.m_values[i];
    }
    return (*this);
}

inline RadMom1D_cState_First_Order& RadMom1D_cState_First_Order::operator *=(const double &b) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] *= b;
  }
  return (*this);
}

inline RadMom1D_cState_First_Order& RadMom1D_cState_First_Order::operator *=(const RadMom1D_cState_First_Order &U) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] *= U.m_values[i];
  }
  return (*this);
}

inline RadMom1D_cState_First_Order& RadMom1D_cState_First_Order::operator /=(const double &b) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] /= b;
  }
  return (*this);
}

inline RadMom1D_cState_First_Order& RadMom1D_cState_First_Order::operator /=(const RadMom1D_cState_First_Order &U) {
  for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++ ) {
    m_values[i] /= U.m_values[i];
  }
  return (*this);
}

/********************************************************
 * RadMom1D_cState_First_Order -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const RadMom1D_cState_First_Order &U) {
    out_file.setf(ios::scientific);

    out_file << " " << U[1] << " " << U[2] << "\n";

    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream &operator >> (istream &in_file, RadMom1D_cState_First_Order &U) {
    in_file.setf(ios::skipws);

    in_file >> U[1] >> U[2];

    in_file.unsetf(ios::skipws);
    return (in_file);
}

/********************************************************
 * RadMom1D_pState_First_Order::RadMom1D_pState_First_Order -- Constructor.       *
 ********************************************************/
inline RadMom1D_pState_First_Order::RadMom1D_pState_First_Order(const RadMom1D_cState_First_Order &U) {
    Nullify();
    Allocate();

    if (U.I0() == ZERO) {
      m_values[0] = U.I0();
      m_values[1] = ZERO;
    } else {
      m_values[0] = U.I0();
      m_values[1] = U.I1x()/U.I0();
    }
}

/********************************************************
 * RadMom1D_pState_First_Order::U -- Conserved solution state.       *  (Start adding closure_type here)
 ********************************************************/
inline double RadMom1D_pState_First_Order :: U(const int &index_U) {
    double U_val;
    switch (index_U) {
        case 1:
            U_val = I0();
            break;
        case 2:
            U_val = I0()*N1x();
            break;
            break;
        default:
            cout << "Incorrect value for index_U = " << index_U << endl;
            exit(0);
            break;
    };
    
    return U_val;
}

inline double RadMom1D_pState_First_Order :: U(const int &index_U) const {
    double U_val;
    switch (index_U) {
        case 1:
            U_val = I0();
            break;
        case 2:
            U_val = I0()*N1x();;
            break;
        default:
            cout << "Incorrect value for index_U = " << index_U << endl;
            exit(0);
            break;
    };
    
    return U_val;
}

inline RadMom1D_cState_First_Order RadMom1D_pState_First_Order::U(void) {
    RadMom1D_cState_First_Order U_temp;

    U_temp[1] = I0();
    U_temp[2] = I0()*N1x();

    return U_temp;
}

inline RadMom1D_cState_First_Order RadMom1D_pState_First_Order::U(void) const {
    RadMom1D_cState_First_Order U_temp;

    U_temp[1] = I0();
    U_temp[2] = I0()*N1x();

    return U_temp;
}

inline RadMom1D_cState_First_Order RadMom1D_pState_First_Order::U(const RadMom1D_pState_First_Order &W) {
    return W.U();
}

inline RadMom1D_cState_First_Order U(const RadMom1D_pState_First_Order &W) {
    return W.U();
}
 
inline void RadMom1D_pState_First_Order :: Compute_Correction_Factor_Roe(double *Correction_Factor,
                                                                         const RadMom1D_pState_First_Order &Wl,
                                                                         const RadMom1D_pState_First_Order &Wr) {
    RadMom1D_cState_First_Order dUrl;
    static double A_Roe[STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER];

    // Compute the jump in the solution between the left and the right states
    dUrl = Wr.U() - Wl.U();

    Compute_Correction_Factor_Approximate_Jacobian_Roe(A_Roe, Wl, Wr);

    for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++) {
      Correction_Factor[i] = ZERO;
      for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; j++) {
        Correction_Factor[i] += A_Roe[i][j] * dUrl.m_values[j];
      }
    }
}
 
inline void RadMom1D_pState_First_Order :: Compute_Correction_Factor_Roe(double *Correction_Factor,
                                                                         const RadMom1D_pState_First_Order &Wl,
                                                                         const RadMom1D_pState_First_Order &Wr) const {
    RadMom1D_cState_First_Order dUrl;
    static double A_Roe[STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER];

    // Compute the jump in the solution between the left and the right states
    dUrl = Wr.U() - Wl.U();

    Compute_Correction_Factor_Approximate_Jacobian_Roe(A_Roe, Wl, Wr);

    for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++) {
      Correction_Factor[i] = ZERO;
      for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; j++) {
        Correction_Factor[i] += A_Roe[i][j] * dUrl.m_values[j];
      }
    }
}
 
inline void RadMom1D_pState_First_Order :: Compute_Correction_Factor_Approximate_Jacobian_Roe(double A_Roe[][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER],
                                                                                              const RadMom1D_pState_First_Order &Wl,
                                                                                              const RadMom1D_pState_First_Order &Wr) {
    Eigenstructure_M1 Eig_M1;
    RadMom1D_pState_First_Order Wstar;
    double lambda_val, lc_val, rc_val;
    static double A_Roe_Temp[STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER];

    Wstar.RoeAverage(Wl, Wr);

    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            Setup_Eigenstructure_P1(Eig_M1);
            break;
        case MOMENT_CLOSURE_M1 :
            // Precompute eigenstructure of the non-gray M1 closure
            if (Wstar.I0() > TOLER_RADIATIVE_DENSITY_M1 &&
                fabs(Wr.I0() - Wl.I0()) > TOLER_RADIATIVE_DENSITY_M1) {
                Setup_Eigenstructure_M1(Eig_M1, Wl, Wr);
            } else {
                Wstar.Setup_Eigenstructure_M1(Eig_M1);
            }
            break;
        default:
            cout << "Closure type not specified" << endl;
            exit(0);
            break;
    }

    for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++) {
        for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; j++) {
            rc_val = Eig_M1.rc_vec[i][j];
            lambda_val = Eig_M1.lambdas[j];
            A_Roe_Temp[i][j] = rc_val*fabs(lambda_val);
        }
    }

    for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++) {
        for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; j++) {
            A_Roe[i][j] = ZERO;
            for (int k = 0; k < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; k++) {
                lc_val = Eig_M1.lc_vec[k][j];
                A_Roe[i][j] += A_Roe_Temp[i][k]*lc_val;
            }
        }
    }
}
 
inline void RadMom1D_pState_First_Order :: Compute_Correction_Factor_Approximate_Jacobian_Roe(double A_Roe[][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER],
                                                                                              const RadMom1D_pState_First_Order &Wl,
                                                                                              const RadMom1D_pState_First_Order &Wr) const {
    Eigenstructure_M1 Eig_M1;
    RadMom1D_pState_First_Order Wstar;
    double lambda_val, lc_val, rc_val;
    static double A_Roe_Temp[STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER][STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER];

    Wstar.RoeAverage(Wl, Wr);

    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            Setup_Eigenstructure_P1(Eig_M1);
            break;
        case MOMENT_CLOSURE_M1 :
            // Precompute eigenstructure of the non-gray M1 closure
            if (Wstar.I0() > TOLER_RADIATIVE_DENSITY_M1 &&
                fabs(Wr.I0() - Wl.I0()) > TOLER_RADIATIVE_DENSITY_M1) {
                Setup_Eigenstructure_M1(Eig_M1, Wl, Wr);
            } else {
                Wstar.Setup_Eigenstructure_M1(Eig_M1);
            }
            break;
        default:
            cout << "Closure type not specified" << endl;
            exit(0);
            break;
    }

    for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++) {
        for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; j++) {
            rc_val = Eig_M1.rc_vec[i][j];
            lambda_val = Eig_M1.lambdas[j];
            A_Roe_Temp[i][j] = rc_val*fabs(lambda_val);
        }
    }

    for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; i++) {
        for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; j++) {
            A_Roe[i][j] = ZERO;
            for (int k = 0; k < STATIC_NUM_VAR_RADMOM1D_FIRST_ORDER; k++) {
                lc_val = Eig_M1.lc_vec[k][j];
                A_Roe[i][j] += A_Roe_Temp[i][k]*lc_val;
            }
        }
    }
}

/********************************************************
 * RadMom1D_pState_First_Order::Fx -- Solution flux (x-direction).   *
 ********************************************************/
inline RadMom1D_cState_First_Order RadMom1D_pState_First_Order::Fx(void) {
    RadMom1D_cState_First_Order flux;

    if (I0() > TOLER_RADIATIVE_DENSITY) {
      flux[1] = I0()*N1x();
      flux[2] = I0()*chi();
    } else {
      flux.Vacuum();
    }
    
    return flux;
}

inline RadMom1D_cState_First_Order RadMom1D_pState_First_Order::Fx(void) const {
    RadMom1D_cState_First_Order flux;

    if (I0() > TOLER_RADIATIVE_DENSITY) {
        flux[1] = I0()*N1x();
        flux[2] = I0()*chi();
    } else {
      flux.Vacuum();
    }
    
    return flux;
}

inline RadMom1D_cState_First_Order RadMom1D_pState_First_Order::Fx(const RadMom1D_pState_First_Order &W) {
    return W.Fx();
}

inline RadMom1D_cState_First_Order Fx(const RadMom1D_pState_First_Order &W) {
    return W.Fx();
}
inline void RadMom1D_pState_First_Order::dFxdU(Eigenstructure_M1 &Eig_M1) {
    double E_dChi2_dE, dChi2_df;

    E_dChi2_dE = Eig_M1.E_dChi2_dE;
    dChi2_df = Eig_M1.dChi2_df;

    Eig_M1.dFdU[0][0] = ZERO;
    Eig_M1.dFdU[0][1] = ONE;
    Eig_M1.dFdU[1][0] = chi() + E_dChi2_dE;
    Eig_M1.dFdU[1][1] = dChi2_df;
}

inline void RadMom1D_pState_First_Order::dFxdU(Eigenstructure_M1 &Eig_M1) const {
    double E_dChi2_dE, dChi2_df;

    E_dChi2_dE = Eig_M1.E_dChi2_dE;
    dChi2_df = Eig_M1.dChi2_df;

    Eig_M1.dFdU[0][0] = ZERO;
    Eig_M1.dFdU[0][1] = ONE;
    Eig_M1.dFdU[1][0] = chi() + E_dChi2_dE;
    Eig_M1.dFdU[1][1] = dChi2_df;
}

/************************************************************
 * RadMom1D_pState_First_Order::lambda_x -- Eigenvalue(s) (x-direction). *
 ************************************************************/
inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order::lambda_x(void) {
    RadMom1D_pState_First_Order lambda;
    
    Evaluate_Lambdas_x(lambda);
    
    return lambda;
}

inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order::lambda_x(void) const {
    RadMom1D_pState_First_Order lambda;
    
    Evaluate_Lambdas_x(lambda);
    
    return lambda;
}

inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order::lambda_x(const RadMom1D_pState_First_Order &W) {
    return W.lambda_x();
}

inline RadMom1D_pState_First_Order lambda_x(const RadMom1D_pState_First_Order &W) {
    return W.lambda_x();
}

//******************************************************************************
// RadMom1D_pState_First_Order::Setup_Eigenstructure_P1.
// This routine sets up the parameters required to compute the eigenvalues and 
// eigenvectors of the flux-Jacobian matrix of the system of moment equations 
// up to first-order resulting from our interpolative-based non-gray M1 closure
//******************************************************************************
inline void RadMom1D_pState_First_Order :: Setup_Eigenstructure_P1(Eigenstructure_M1 &Eig_P1) {
    // Set entries of the arrays contained in the structure Eig_M1 to zero
    Eig_P1.Nullify();

    // Compute the flux-Jacobian of the resulting closed system of moment equations
    dFxdU(Eig_P1);
    
    // Compute the eigenvalues of the the flux-Jacobian matrix
    Evaluate_Lambdas_x(Eig_P1);
    
    // Compute the eigenvectors of the the flux-Jacobian matrix
    // Right eigenvectors
    rc(Eig_P1);
    // Left eigenvectors
    lc(Eig_P1);
}

inline void RadMom1D_pState_First_Order :: Setup_Eigenstructure_P1(Eigenstructure_M1 &Eig_P1) const {
    // Set entries of the arrays contained in the structure Eig_M1 to zero
    Eig_P1.Nullify();

    // Compute the flux-Jacobian of the resulting closed system of moment equations
    dFxdU(Eig_P1);
    
    // Compute the eigenvalues of the the flux-Jacobian matrix
    Evaluate_Lambdas_x(Eig_P1);
    
    // Compute the eigenvectors of the the flux-Jacobian matrix
    // Right eigenvectors
    rc(Eig_P1);
    // Left eigenvectors
    lc(Eig_P1);
}

//******************************************************************************
// RadMom1D_pState_First_Order::Setup_Eigenstructure_M1.
// This routine sets up the parameters required to compute the eigenvalues and 
// eigenvectors of the flux-Jacobian matrix of the system of moment equations 
// up to first-order resulting from our interpolative-based non-gray M1 closure
//******************************************************************************
inline void RadMom1D_pState_First_Order :: Setup_Eigenstructure_M1(Eigenstructure_M1 &Eig_M1) {
    double E_dChi2_dE, dChi2_df;

    // Set entries of the arrays contained in the structure Eig_M1 to zero
    Eig_M1.Nullify();
    
    // Compute derivatives of the second-order closing fluxes with respect to the
    // lower-order angular moments
    Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_df);
    
    // Store the value in the Eigenstructure_M1 data structure
    Eig_M1.E_dChi2_dE = E_dChi2_dE;
    Eig_M1.dChi2_df = dChi2_df;
    
    // Compute the flux-Jacobian of the resulting closed system of moment equations
    dFxdU(Eig_M1);
    
    // Compute the eigenvalues of the the flux-Jacobian matrix
    Evaluate_Lambdas_x(Eig_M1);
    
    // Compute the eigenvectors of the the flux-Jacobian matrix
    // Right eigenvectors
    rc(Eig_M1);
    // Left eigenvectors
    lc(Eig_M1);
}

inline void RadMom1D_pState_First_Order :: Setup_Eigenstructure_M1(Eigenstructure_M1 &Eig_M1) const {
    double E_dChi2_dE, dChi2_df;
    
    // Set entries of the arrays contained in the structure Eig_M1 to zero
    Eig_M1.Nullify();

    // Compute derivatives of the second-order closing fluxes with respect to the
    // lower-order angular moments
    Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_df);
    
    // Store the value in the Eigenstructure_M1 data structure
    Eig_M1.E_dChi2_dE = E_dChi2_dE;
    Eig_M1.dChi2_df = dChi2_df;
    
    // Compute the flux-Jacobian of the resulting closed system of moment equations
    dFxdU(Eig_M1);
    
    // Compute the eigenvalues of the the flux-Jacobian matrix
    Evaluate_Lambdas_x(Eig_M1);
    
    // Compute the eigenvectors of the the flux-Jacobian matrix
    // Right eigenvectors
    rc(Eig_M1);
    // Left eigenvectors
    lc(Eig_M1);
}

inline void RadMom1D_pState_First_Order :: Setup_Eigenstructure_M1(Eigenstructure_M1 &Eig_M1, const RadMom1D_pState_First_Order &Wl, const RadMom1D_pState_First_Order &Wr) {
    RadMom1D_pState_First_Order Wl_Roe, Wr_Roe, Wstar_Roe, Wstar;
    RadMom1D_cState_First_Order Ustar;
    RadMom1D_cState_First_Order Ul, Ur;
    
    Ul = Wl.U();
    Ur = Wr.U();
    
    Wl_Roe = U_to_W_Roe(Ul);
    Wr_Roe = U_to_W_Roe(Ur);
    Wstar_Roe.AverageStates(Wl_Roe, Wr_Roe);
    
    Ustar = W_Roe_to_U(Wstar_Roe);
    Wstar = Ustar.W();
    
    // Set entries of the arrays contained in the structure Eig_M1 to zero
    Eig_M1.Nullify();
    
    // Compute the Roe matrix
    Approximate_Roe_Matrix(Eig_M1, Wl, Wr);
    
    // Compute the eigenvalues of the the flux-Jacobian matrix
    Wstar.Evaluate_Lambdas_x(Eig_M1);
    
    // Compute the eigenvectors of the the flux-Jacobian matrix
    // Right eigenvectors
    Wstar.rc(Eig_M1);
    // Left eigenvectors
    Wstar.lc(Eig_M1);
}

inline void RadMom1D_pState_First_Order :: Setup_Eigenstructure_M1(Eigenstructure_M1 &Eig_M1, const RadMom1D_pState_First_Order &Wl, const RadMom1D_pState_First_Order &Wr) const {
    RadMom1D_pState_First_Order Wl_Roe, Wr_Roe, Wstar_Roe, Wstar;
    RadMom1D_cState_First_Order Ustar;
    RadMom1D_cState_First_Order Ul, Ur;
    
    Ul = Wl.U();
    Ur = Wr.U();
    
    Wl_Roe = U_to_W_Roe(Ul);
    Wr_Roe = U_to_W_Roe(Ur);
    Wstar_Roe.AverageStates(Wl_Roe, Wr_Roe);
    
    Ustar = W_Roe_to_U(Wstar_Roe);
    Wstar = Ustar.W();
    
    // Set entries of the arrays contained in the structure Eig_M1 to zero
    Eig_M1.Nullify();
    
    // Compute the Roe matrix
    Approximate_Roe_Matrix(Eig_M1, Wl, Wr);
    
    // Compute the eigenvalues of the the flux-Jacobian matrix
    Wstar.Evaluate_Lambdas_x(Eig_M1);
    
    // Compute the eigenvectors of the the flux-Jacobian matrix
    // Right eigenvectors
    Wstar.rc(Eig_M1);
    // Left eigenvectors
    Wstar.lc(Eig_M1);
}

inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order :: U_to_W_Roe(const RadMom1D_cState_First_Order &U) {
    RadMom1D_pState_First_Order W_Roe;

    W_Roe[1] = sqrt(U.I0());
    W_Roe[2] = U.I1x()/sqrt(U.I0());

    return W_Roe;
}

inline RadMom1D_pState_First_Order RadMom1D_pState_First_Order :: U_to_W_Roe(const RadMom1D_cState_First_Order &U) const {
    RadMom1D_pState_First_Order W_Roe;

    W_Roe[1] = sqrt(U.I0());
    W_Roe[2] = U.I1x()/sqrt(U.I0());

    return W_Roe;
}

inline RadMom1D_cState_First_Order RadMom1D_pState_First_Order :: W_Roe_to_U(const RadMom1D_pState_First_Order &W_Roe) {
    RadMom1D_cState_First_Order U;

    U[1] = pow(W_Roe.I0(), 2);
    U[2] = W_Roe.I0()*W_Roe.N1x();

    return U;
}

inline RadMom1D_cState_First_Order RadMom1D_pState_First_Order :: W_Roe_to_U(const RadMom1D_pState_First_Order &W_Roe) const {
    RadMom1D_cState_First_Order U;

    U[1] = pow(W_Roe.I0(), 2);
    U[2] = W_Roe.I0()*W_Roe.N1x();

    return U;
}

inline void RadMom1D_pState_First_Order :: Approximate_Roe_Matrix(Eigenstructure_M1 &Eig_M1,
                                                           const RadMom1D_pState_First_Order &Wl,
                                                           const RadMom1D_pState_First_Order &Wr) {
    double chi2_l, chi2_r;
    double C, D;
    chi2_l = Wl.chi();
    chi2_r = Wr.chi();
    C = (chi2_r - chi2_l)/(Wr.N1x() - Wl.N1x());
    D = (chi2_l*Wr.N1x() - chi2_r*Wl.N1x())/(Wr.N1x() - Wl.N1x());
    Eig_M1.dFdU[0][0] = ZERO;
    Eig_M1.dFdU[0][1] = ONE;
    Eig_M1.dFdU[1][0] = D;
    Eig_M1.dFdU[1][1] = C;
}

inline void RadMom1D_pState_First_Order :: Approximate_Roe_Matrix(Eigenstructure_M1 &Eig_M1,
                                                           const RadMom1D_pState_First_Order &Wl,
                                                           const RadMom1D_pState_First_Order &Wr) const {
    double chi2_l, chi2_r;
    double C, D;
    chi2_l = Wl.chi();
    chi2_r = Wr.chi();
    C = (chi2_r - chi2_l)/(Wr.N1x() - Wl.N1x());
    D = (chi2_l*Wr.N1x() - chi2_r*Wl.N1x())/(Wr.N1x() - Wl.N1x());
    Eig_M1.dFdU[0][0] = ZERO;
    Eig_M1.dFdU[0][1] = ONE;
    Eig_M1.dFdU[1][0] = D;
    Eig_M1.dFdU[1][1] = C;
}

//******************************************************************************
//RadMom1D_pState_First_Order::rc -- Conserved right eigenvector (x-direction).                 
//******************************************************************************
inline void RadMom1D_pState_First_Order::rc(Eigenstructure_M1 &Eig_M1) {
    // Find the vector rc such that: A rc = \lam_i rc
    // where A is the flux Jacobian and \lam_i is the an eigenvalue of the matrix A
    // Right eigenvectors are stored columnwise in the matrix of eigenvectors

    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            // First column
            Eig_M1.rc_vec[0][0] = -sqrt(THREE);
            Eig_M1.rc_vec[1][0] = ONE;

            // Second column
            Eig_M1.rc_vec[0][1] = sqrt(THREE);
            Eig_M1.rc_vec[1][1] = ONE;
            break;
        case MOMENT_CLOSURE_M1:
            // First column
            Eig_M1.rc_vec[0][0] = ONE;
            Eig_M1.rc_vec[1][0] = Eig_M1.lambdas[0];

            // Second column
            Eig_M1.rc_vec[0][1] = ONE;
            Eig_M1.rc_vec[1][1] = Eig_M1.lambdas[1];
            break;
        default :
            cout << "closure_type not specified !!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    };
    // (A - lam I) v = 0
    // [    -lam                        1    ] [v1]    [0]
    //                                           =
    // [E_dChi2_dE - N1x*dChi2df      dChi2df] [v2]    [0]
    // -lam*v1 + v2 = 0
    // (E_dChi2_dE - N1x*dChi2df)*v1 + dChi2df*v2 = 0

    // v2 = lam*v1
    // v2 = - (E_dChi2_dE - N1x*dChi2df)*v1/dChi2_df
    //      [ 1         1 ]
    // rc =
    //      [lam1     lam2]
}

inline void RadMom1D_pState_First_Order::rc(Eigenstructure_M1 &Eig_M1) const {
    // Find the vector rc such that: A rc = \lam_i rc
    // where A is the flux Jacobian and \lam_i is the an eigenvalue of the matrix A
    // Right eigenvectors are stored columnwise in the matrix of eigenvectors

    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            // First column
            Eig_M1.rc_vec[0][0] = -sqrt(THREE);
            Eig_M1.rc_vec[1][0] = ONE;

            // Second column
            Eig_M1.rc_vec[0][1] = sqrt(THREE);
            Eig_M1.rc_vec[1][1] = ONE;
            break;
        case MOMENT_CLOSURE_M1:
            // First column
            Eig_M1.rc_vec[0][0] = ONE;
            Eig_M1.rc_vec[1][0] = Eig_M1.lambdas[0];

            // Second column
            Eig_M1.rc_vec[0][1] = ONE;
            Eig_M1.rc_vec[1][1] = Eig_M1.lambdas[1];
            break;
        default :
            cout << "closure_type not specified !!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    };
    // (A - lam I) v = 0
    // [    -lam                        1    ] [v1]    [0]
    //                                           =
    // [E_dChi2_dE - N1x*dChi2df      dChi2df] [v2]    [0]
    // -lam*v1 + v2 = 0
    // (E_dChi2_dE - N1x*dChi2df)*v1 + dChi2df*v2 = 0

    // v2 = lam*v1
    // v2 = - (E_dChi2_dE - N1x*dChi2df)*v1/dChi2_df
    //      [ 1         1 ]
    // rc =
    //      [lam1     lam2]
}

inline void rc(const RadMom1D_pState_First_Order &W, Eigenstructure_M1 &Eig_M1) {
    W.rc(Eig_M1);
}

//******************************************************************************
//RadMom1D_pState_First_Order::rc -- Conserved left eigenvector (x-direction).                 
//******************************************************************************
inline void RadMom1D_pState_First_Order::lc(Eigenstructure_M1 &Eig_M1) {
    double deter;
    // Find the vector rc such that: A rc = \lam_i rc
    // where A is the flux Jacobian and \lam_i is the an eigenvalue of the matrix A
    // Left eigenvectors are stored row-wise in the matrix of eigenvectors

    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            // First row
            Eig_M1.lc_vec[0][0] = -sqrt(THREE)/SIX;
            Eig_M1.lc_vec[0][1] = HALF;

            // Second row
            Eig_M1.lc_vec[1][0] = sqrt(THREE)/SIX;
            Eig_M1.lc_vec[1][1] = HALF;
            break;
        case MOMENT_CLOSURE_M1:
            // First row
            deter = Eig_M1.lambdas[1] - Eig_M1.lambdas[0];
            Eig_M1.lc_vec[0][0] = Eig_M1.lambdas[1]/deter;
            Eig_M1.lc_vec[1][0] = -Eig_M1.lambdas[0]/deter;

            // Second row
            Eig_M1.lc_vec[0][1] = -ONE/deter;
            Eig_M1.lc_vec[1][1] = ONE/deter;
            break;
        default :
            cout << "closure_type not specified !!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    };

    // lc (A - lam I) = 0
    //            [    -lam                        1    ]    [0]
    // [lc1  lc2]                                          =
    //            [E_dChi2_dE - N1x*dChi2df      dChi2df]    [0]
    // -lam*lc1 + (E_dChi2_dE - N1x*dChi2df)*lc2 = 0
    // lc1 + dChi2df*lc2 = 0

    // lc2 = lam*lc1/(E_dChi2_dE - N1x*dChi2df);
    // lc2 = -lc1/dChi2_df
    //      [ 1         1 ]
    // rc =
    //      [lam1     lam2]
}

inline void RadMom1D_pState_First_Order::lc(Eigenstructure_M1 &Eig_M1) const {
    double deter;
    // Find the vector rc such that: A rc = \lam_i rc
    // where A is the flux Jacobian and \lam_i is the an eigenvalue of the matrix A
    // Left eigenvectors are stored row-wise in the matrix of eigenvectors

    switch(closure_type) {
        case MOMENT_CLOSURE_P1 :
            // First row
            Eig_M1.lc_vec[0][0] = -sqrt(THREE)/SIX;
            Eig_M1.lc_vec[0][1] = HALF;

            // Second row
            Eig_M1.lc_vec[1][0] = sqrt(THREE)/SIX;
            Eig_M1.lc_vec[1][1] = HALF;
            break;
        case MOMENT_CLOSURE_M1:
            // First row
            deter = Eig_M1.lambdas[1] - Eig_M1.lambdas[0];
            Eig_M1.lc_vec[0][0] = Eig_M1.lambdas[1]/deter;
            Eig_M1.lc_vec[1][0] = -Eig_M1.lambdas[0]/deter;

            // Second row
            Eig_M1.lc_vec[0][1] = -ONE/deter;
            Eig_M1.lc_vec[1][1] = ONE/deter;
            break;
        default :
            cout << "closure_type not specified !!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    };
}

inline void lc(const RadMom1D_pState_First_Order &W, Eigenstructure_M1 &Eig_M1) {
    W.lc(Eig_M1);
}

/**********************************************************************
 * RadMom1D_pState_First_Order::S -- Include all source term vectors and    *
 *                             Jacobians.                             *
 * Regular Source Term
**********************************************************************/
inline RadMom1D_cState_First_Order  RadMom1D_cState_First_Order::S(const Medium1D_State &M ) {
    RadMom1D_cState_First_Order Source;

    Source[1] = M.kappa()*(FOUR*PI*M.Ib() - I0());
    Source[2] = -M.beta()*I1x();
  
  return Source;
}

inline double RadMom1D_cState_First_Order::Sr(const Medium1D_State &M ) {
    double Source_term;
    Source_term = M.kappa()*(FOUR*PI*M.Ib() - I0());
  
    return -Source_term;
}

inline double RadMom1D_cState_First_Order::Sr(const Medium1D_State &M ) const {
    double Source_term;
    Source_term = M.kappa()*(FOUR*PI*M.Ib() - I0());

    return -Source_term;
}

/*****************************************************************
 * First order piecewise linear solution reconstruction
 ******************************************************************/
inline void RadMom1D_pState_First_Order::Reconstruct( const RadMom1D_pState_First_Order &Wc, 
                                                      const RadMom1D_pState_First_Order &phi, 
                                                      const RadMom1D_pState_First_Order &dWdx,
                                                      const double &dX) {
    m_values[0] = Wc.I0() + phi.I0()*dWdx.I0()*dX;
    m_values[1] = Wc.N1x() + phi.N1x()*dWdx.N1x()*dX;
}

/********************************************************
 * RadMom1D_cState_First_Order::RadMom1D_cState_First_Order -- Constructor.     *
 ********************************************************/
inline RadMom1D_cState_First_Order::RadMom1D_cState_First_Order(const RadMom1D_pState_First_Order &W) {
    Nullify();
    Allocate();
    m_values[0] = W.I0(); 
    m_values[1] = W.I0()*W.N1x();
}

/********************************************************
 * RadMom1D_cState_First_Order::W -- Primitive solution state.       *
 ********************************************************/
inline RadMom1D_pState_First_Order RadMom1D_cState_First_Order::W(void) {
    RadMom1D_pState_First_Order W_temp;
    
    W_temp[1] = I0();
    W_temp[2] = I1x()/I0();
    
    return W_temp;
}

inline RadMom1D_pState_First_Order RadMom1D_cState_First_Order::W(void) const {
    RadMom1D_pState_First_Order W_temp;
    
    W_temp[1] = I0();
    W_temp[2] = I1x()/I0();

    return W_temp;
}

inline RadMom1D_pState_First_Order RadMom1D_cState_First_Order::W(const RadMom1D_cState_First_Order &U) {
    return U.W();
}

inline RadMom1D_pState_First_Order W(const RadMom1D_cState_First_Order &U) {
    return U.W();
}

/********************************************************
 * RadMom1D_cState_First_Order::Fx -- Solution flux (x-direction).   *
 ********************************************************/
inline RadMom1D_cState_First_Order RadMom1D_cState_First_Order::Fx(void) {
    return (W().Fx());
}

inline RadMom1D_cState_First_Order RadMom1D_cState_First_Order::Fx(void) const {
    return (W().Fx());
}

inline RadMom1D_cState_First_Order RadMom1D_cState_First_Order::Fx(const RadMom1D_cState_First_Order &U) {
    return (U.W().Fx());
}

inline RadMom1D_cState_First_Order Fx(const RadMom1D_cState_First_Order &U) {
    return (U.W().Fx());
}

#endif /* _RADMOM1D_STATE_FIRST_ORDER_INCLUDED  */
