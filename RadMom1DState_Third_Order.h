/*!\file RadMom1DState_Third_Order.h
  \brief Header file defining 1D RadMom Solution State Classes. */

#ifndef _RADMOM1D_STATE_THIRD_ORDER_INCLUDED
#define _RADMOM1D_STATE_THIRD_ORDER_INCLUDED

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

#ifndef _RADMOM1D_STATE_INCLUDED
#include "RadMom1DState.h"
#endif //_RADMOM1D_STATE_INCLUDED

/* Define some constants. */
#define STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER 4

/* Define additional structures. */
struct Eigenstructure_P3 {
    double lambdas[STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER];
    double rc_vec[STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER][STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER];
    double lc_vec[STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER][STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER];
};

/* Define the classes. */
class RadMom1D_cState_Third_Order;

/*!
 * Class: RadMom1D_cState_Third_Order
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
 *     a        -- Return 8*pi^5*k^4/(15*h^3*c^3).
 *     U        -- Return conserved solution state.
 *     Fx       -- Return x-direction solution flux.
 *     dFxdU    -- Return x-direction flux Jacobian.
 *     lambda_x -- Return x-direction eigenvalue(s).
 *     rp       -- Return primitive right eigenvector (x-direction).
 *     rc       -- Return conserved right eigenvector (x-direction).
 *     S        -- Return source term vector.
 *     Sr        -- Return radiative source term.
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
class RadMom1D_pState_Third_Order : public RadMom1D_pState<RadMom1D_cState_Third_Order, RadMom1D_pState_Third_Order>{
protected:
  
private:
public:
  #ifdef STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER
  double  m_values[STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER];
  #else 
  double* m_values;
  #endif
  
    //@{ @name constants:
    static int    NUM_VAR_RADMOM1D_THIRD_ORDER; //!< total number of RadMom1D state variables
    //@}
    
    //!

  void Setup_Eigenstructure_P3(Eigenstructure_P3 &Eig_P3) const;
  
  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  RadMom1D_pState_Third_Order(void) {
      Nullify();
      Allocate();
      
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
          m_values[i] = ZERO;
      }
  }

  //! Copy constructor.
  RadMom1D_pState_Third_Order(const RadMom1D_pState_Third_Order &W) {
      Nullify();
      Allocate();
      
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
          m_values[i] = W.m_values[i];  
      }
  }

  //! Assignment constructor.
  RadMom1D_pState_Third_Order(const RadMom1D_cState_Third_Order &U);

  //! Assignment constructor.
  RadMom1D_pState_Third_Order(const double &den, const double &fx,
                              const double &p_xx, const double &q_xxx) {
      Nullify();
      Allocate();
      
      m_values[0] = den;
      m_values[1] = fx; 
      m_values[2] = p_xx; 
      m_values[3] = q_xxx; 
  }
  
  //! Value Constructor
  explicit RadMom1D_pState_Third_Order(const double &Val);

  /* Destructor. */
  ~RadMom1D_pState_Third_Order(void) { Deallocate(); }
  // Use automatically generated destructor.
  //@}

  //@{ @name Useful operators.
  //! Return the number of variables.
  
  static int StaticNumVar(void) { return STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; }
  static int NumVar(void) { return NUM_VAR_RADMOM1D_THIRD_ORDER; }
  
  //! Copy operator.
  void Copy(const RadMom1D_pState_Third_Order &W) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
          m_values[i] = W.m_values[i];
      }
  }
  
  void Copy_to_W(RadMom1D_pState_Third_Order &W) const {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
          W.m_values[i] = m_values[i];
      }
  }
  
  //! Set array pointers to null
  void Nullify();

  //! Vacuum operator.
  void Vacuum(void) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
          m_values[i] = ZERO;
      }
  }

  //! One operator. Set the solution to ONE.
  void Ones(void) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
          m_values[i] = ONE;
      }
  }
  
  double I0() const {return m_values[0];}
  double N1x() const {return m_values[1];}
  double N2xx() const {return m_values[2];}
  double N3xxx() const {return m_values[3];}
  double r_xxxx() const {return ((SIX/SEVEN)*N2xx()-(THREE/THIRTY_FIVE));}
  
  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const {
      return 0;
  }
  //@}

  //@{ @name State functions.
  //! Absolute flux (scalar).
  double fsca() const;

 //! Radiative flux squared (scalar).
  double f2() const;
  
  //@{ @name Conserved solution state.
  double U(const int &index_U) const;
  RadMom1D_cState_Third_Order U(void) const;
  RadMom1D_cState_Third_Order U(const RadMom1D_pState_Third_Order &W);
  friend RadMom1D_cState_Third_Order U(const RadMom1D_pState_Third_Order &W);
  //@}

  //@{ @name Solution flux and Jacobian (x-direction).
  RadMom1D_cState_Third_Order Fx(void) const;
  RadMom1D_cState_Third_Order Fx(const RadMom1D_pState_Third_Order &W);
  friend RadMom1D_cState_Third_Order Fx(const RadMom1D_pState_Third_Order &W);
  //@}

  //@{ @name Eigenvalue(s) (x-direction).
  RadMom1D_pState_Third_Order lambda_x(void) const;
  RadMom1D_pState_Third_Order lambda_x(const RadMom1D_pState_Third_Order &W);
  friend RadMom1D_pState_Third_Order lambda_x(const RadMom1D_pState_Third_Order &W);
  friend double lambda_x(const RadMom1D_pState_Third_Order &W, int index);
  RadMom1D_pState_Third_Order lambda_x(const double &V) const;
  //@}

  void Compute_Correction_Factor_Roe(double *Correction_Factor,
                                     const RadMom1D_pState_Third_Order &Wl,
                                     const RadMom1D_pState_Third_Order &Wr) const;

  void Compute_Correction_Factor_Approximate_Jacobian_Roe(double A_Roe[][STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER]) const;
                                                          
  //@{ @name Conserved right eigenvector (x-direction).
  double rc(const int &index_U, const int &j) const;
  void rc(Eigenstructure_P3 &Eig_P3) const;
  friend void rc(const RadMom1D_pState_Third_Order &W, Eigenstructure_P3 &Eig_P3);
  //@}

  //@{ @name Primitive right eigenvector (x-direction).
  RadMom1D_pState_Third_Order rp(const int &index_U) const;
  friend RadMom1D_pState_Third_Order rp(const RadMom1D_pState_Third_Order &W, const int &index_U);
  //@}
  
  //@{ @name Conserved left eigenvector (x-direction).
  double lc(const int &index_U, const int &j) const;
  void lc(Eigenstructure_P3 &Eig_P3) const;
  friend void lc(const RadMom1D_pState_Third_Order &W, Eigenstructure_P3 &Eig_P3);
  //@}
  
  //@{ @name Primitive left eigenvector (x-direction).
  RadMom1D_pState_Third_Order lp(const int &index_U) const;
  friend RadMom1D_pState_Third_Order lp(const RadMom1D_pState_Third_Order &W, const int &index_U);
  //@}

  //@{ @name First order linear solution reconstruction
  void Reconstruct( const RadMom1D_pState_Third_Order &Wc, 
                    const RadMom1D_pState_Third_Order &phi, 
                    const RadMom1D_pState_Third_Order &dWdx, 
                    const double &dX);
  
  //! memory allocation / deallocation for the I array
  void Allocate();
  void Deallocate();
  
  //! Setup state static variables
  //! memory allocation / deallocation for the static arrays
  static void SetupStatic( const int &ScatteringFunc,
                           const int &i_Moment_Closure,
                           const int &i_AbsorptionModel);
  static void DeallocateStatic();
  
  void RoeAverage(const RadMom1D_pState_Third_Order &Wl,
                  const RadMom1D_pState_Third_Order &Wr);
  void Rotate(const double &norm_dir);
  void Reflect(RadMom1D_pState_Third_Order W_inner,
               const double &norm_dir);
  
  void Characteristic(RadMom1D_pState_Third_Order W_inner,
                     const double &wall_temperature, 
                     const double &wall_emissivity, 
                     const double &norm_dir);
  void PartialFlux_n(RadMom1D_pState_Third_Order W_inner,
                     const double &wall_temperature, 
                     const double &wall_emissivity, 
                     const double &norm_dir);
  void Gray_Wall(RadMom1D_pState_Third_Order W_wall,
                 const double &wall_temperature, 
                 const double &wall_emissivity, 
                 const double &norm_dir);
  void Marshak_n(RadMom1D_pState_Third_Order W_inner,
                 const double &wall_temperature, 
                 const double &wall_emissivity, 
                 const double &norm_dir);
  void PartialMoments_n(RadMom1D_pState_Third_Order W_inner,
                        const double &wall_temperature, 
                        const double &wall_emissivity, 
                        const double &norm_dir);

  //! Index operator.
  double &operator[](int index) {
      return W_index(*this, index);
  }
  const double &operator[](int index) const {
      return W_index(*this, index);
  }
  
  //@{ @name Binary arithmetic operators.
  RadMom1D_pState_Third_Order operator +(const RadMom1D_pState_Third_Order &W);
  RadMom1D_pState_Third_Order operator -(const RadMom1D_pState_Third_Order &W);
  RadMom1D_pState_Third_Order operator *(const double &b);
  double operator *(const RadMom1D_pState_Third_Order &W);
  friend RadMom1D_pState_Third_Order operator *(const double &b, const RadMom1D_pState_Third_Order &W);
  RadMom1D_pState_Third_Order operator /(const double &b);
  RadMom1D_pState_Third_Order operator /(const RadMom1D_pState_Third_Order &W);
  RadMom1D_pState_Third_Order operator ^(const RadMom1D_pState_Third_Order &W);
  
  /* Assignment operator. */
  RadMom1D_pState_Third_Order &operator = (const RadMom1D_pState_Third_Order &W);
  
  //@{ @name Shortcut arithmetic operators.
  RadMom1D_pState_Third_Order &operator +=(const RadMom1D_pState_Third_Order &W);
  RadMom1D_pState_Third_Order &operator -=(const RadMom1D_pState_Third_Order &W);
  RadMom1D_pState_Third_Order &operator *=(const double &b);
  RadMom1D_pState_Third_Order &operator *=(const RadMom1D_pState_Third_Order &W);
  RadMom1D_pState_Third_Order &operator /=(const double &b);
  RadMom1D_pState_Third_Order &operator /=(const RadMom1D_pState_Third_Order &W);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const RadMom1D_pState_Third_Order &W);
  friend istream &operator >> (istream &in_file,  RadMom1D_pState_Third_Order &W);
  //@}
};

 // end of RadMom1D_pState_Third_Order class

/*!
 * Class: RadMom1D_cState_Third_Order
 *
 * @brief Conserved variable solution state class definition for radiative
 *        transfer.
 *
 * Conserved variable solution state class definition for radiative
 * transfer.
 *     E        -- Return radiative energy (density).
 *     F        -- Return radiative flux (vector).
 *     fsca     -- Return absolute normalized flux (scalar).
 *     f2       -- Return square of radiative flux (scalar).
 *     xi       -- Return sqrt(4-3*fsca^2).
 *     chi      -- Return (3+4*fsca^2)/(5+2*xi).
 *     k        -- Return Boltzmann constant.
 *     h        -- Return Planck's constant.
 *     c        -- Return speed of light.
 *     a        -- Return 8*pi^5*k^4/(15*h^3*c^3).
 *     W        -- Return primitive solution state.
 *     Flux     -- Return x-direction solution flux.
 *     Fx       -- Return x-direction solution flux.
 *     dFxdU    -- Return x-direction flux Jacobian.
 *     S        -- Return source term vector.
 *     Sr        -- Return radiative source term.
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
class RadMom1D_cState_Third_Order : public RadMom1D_cState<class RadMom1D_cState_Third_Order, class RadMom1D_pState_Third_Order> {
protected:
private:
public:
  #ifdef STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER
  double  m_values[STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER];
  #else 
  double* m_values;
  #endif
  
    //@{ @name constants:
  static int    NUM_VAR_RADMOM1D_THIRD_ORDER; //!< total number of RadMom1D state variables
  //@}
  
  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  RadMom1D_cState_Third_Order(void) {
      Nullify();
      Allocate();
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
          m_values[i] = ZERO;
      }
  }

  //! Copy constructor.
  RadMom1D_cState_Third_Order(const RadMom1D_cState_Third_Order &U) {
      Nullify();
      Allocate();
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
          m_values[i] = U.m_values[i];  
      }
  }

  //! Copy constructor.
  RadMom1D_cState_Third_Order(const RadMom1D_pState_Third_Order &W);

  //! Value Constructor
  explicit RadMom1D_cState_Third_Order(const double &Val);

  //! Assignment constructor.
  RadMom1D_cState_Third_Order(const double &den, const double &Fx,
                              const double &P_xx, const double &Q_xxx) {
      Nullify();
      Allocate();
      
      m_values[0] = den;
      m_values[1] = Fx;
      m_values[2] = P_xx; 
      m_values[3] = Q_xxx; 
  }

  /* Destructor. */
  ~RadMom1D_cState_Third_Order(void) { Deallocate(); }
  // Use automatically generated destructor.
  //@}

  //@{ @name Useful operators.
  //! Return the number of variables.
  static int StaticNumVar(void) { return STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; }
  static int NumVar(void) { return NUM_VAR_RADMOM1D_THIRD_ORDER; }

  //! Copy operator.
  void Copy(const RadMom1D_cState_Third_Order &U) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
          m_values[i] = U.m_values[i];
      }
  }
  
  //! Set array pointers to null
  void Nullify();

  //! Vacuum operator.
  void Vacuum(void) {
      for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
          m_values[i] = ZERO;
      }
  }
  
  double I0() const {return m_values[0];}
  double I1x() const {return m_values[1];}
  double I2xx() const {return m_values[2];}
  double I3xxx() const {return m_values[3];}
  double R_xxxx() const {return ((SIX/SEVEN)*I2xx()-(THREE/THIRTY_FIVE)*I0());}
  
  //! Check for unphysical state properties.
  int Unphysical_Properties() const {
      
      return 0;
  }

  //@{ @name State functions.
  //! Absolute radiative flux (scalar).
  double fsca();
  double fsca() const;

  //! Square of radiative flux (scalar).
  double f2();
  double f2() const;
  
    //@{ @name Primitive solution state.
  RadMom1D_pState_Third_Order W(void) const;
  RadMom1D_pState_Third_Order W(const RadMom1D_cState_Third_Order &U);
  friend RadMom1D_pState_Third_Order W(const RadMom1D_cState_Third_Order &U);
  //@}

  //@{ @name Solution flux and Jacobian (x-direction).
  RadMom1D_cState_Third_Order Fx(void) const;
  RadMom1D_cState_Third_Order Fx(const RadMom1D_cState_Third_Order &U);
  friend RadMom1D_cState_Third_Order Fx(const RadMom1D_cState_Third_Order &U);
  //@}

  //@{ @name Include all source vectors and Jacobians.
  RadMom1D_cState_Third_Order S(const Medium1D_State &M );
  double Sr(const Medium1D_State &M);
  
  //! memory allocation / deallocation for the I array
  void Allocate();
  void Deallocate();
  
  //! Setup state static variables
  //! memory allocation / deallocation for the static arrays
  static void SetupStatic( const int &ScatteringFunc,
                           const int &i_Moment_Closure,
                           const int &i_AbsorptionModel);
  static void DeallocateStatic();
  
  void Rotate(const double &norm_dir);
  void Set_ICs(const double &Medium_Temperature);
  void Set_ICs_Beam(const double &Medium_Temperature);
  void Set_ICs_Intensity(const double &Ib_wall);
  void Set_BCs(const double *Intensity, const double & norm_dir);
  
  void Partial_Numerical_Flux(RadMom1D_pState_Third_Order W_inner,
                              const double &norm_dir);
  
  //! Index operator.
  double &operator[](int index) {
      return U_index(*this, index);
  }
  const double &operator[](int index) const {
    return U_index(*this, index);
  }

  //@{ @name Binary arithmetic operators.
  RadMom1D_cState_Third_Order operator +(const RadMom1D_cState_Third_Order &U);
  RadMom1D_cState_Third_Order operator -(const RadMom1D_cState_Third_Order &U);
  RadMom1D_cState_Third_Order operator *(const double &b);
  double operator *(const RadMom1D_cState_Third_Order &U);
  friend RadMom1D_cState_Third_Order operator *(const double &b, const RadMom1D_cState_Third_Order &U);
  RadMom1D_cState_Third_Order operator /(const double &b);
  RadMom1D_cState_Third_Order operator /(const RadMom1D_cState_Third_Order &U);
  RadMom1D_cState_Third_Order operator ^(const RadMom1D_cState_Third_Order &U);
  
  /* Assignment operator. */
  RadMom1D_cState_Third_Order &operator = (const RadMom1D_cState_Third_Order &U);
  
  //@{ @name Shortcut arithmetic operators.
  RadMom1D_cState_Third_Order &operator +=(const RadMom1D_cState_Third_Order &U);
  RadMom1D_cState_Third_Order &operator -=(const RadMom1D_cState_Third_Order &U);
  RadMom1D_cState_Third_Order &operator *=(const double &b);
  RadMom1D_cState_Third_Order &operator *=(const RadMom1D_cState_Third_Order &U);
  RadMom1D_cState_Third_Order &operator /=(const double &b);
  RadMom1D_cState_Third_Order &operator /=(const RadMom1D_cState_Third_Order &U);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const RadMom1D_cState_Third_Order &U);
  friend istream &operator >> (istream &in_file,  RadMom1D_cState_Third_Order &U);
  //@}

};

// End of RadMom1D_cState_Third_Order class

/****************************************************************//**
 * RadMom1D_pState_Third_Order -- Binary arithmetic operators.
 ********************************************************************/
//! addition
inline RadMom1D_pState_Third_Order RadMom1D_pState_Third_Order :: operator +(const RadMom1D_pState_Third_Order &W) {
  RadMom1D_pState_Third_Order W_temp(*this);
  addition(W_temp, W);
  return W_temp;
}

//! subtraction
inline RadMom1D_pState_Third_Order RadMom1D_pState_Third_Order :: operator -(const RadMom1D_pState_Third_Order &W) {
  RadMom1D_pState_Third_Order W_temp(*this);
  subtraction(W_temp, W);
  return W_temp;
}

//! scalar multiplication
inline RadMom1D_pState_Third_Order RadMom1D_pState_Third_Order :: operator *(const double &b) {
  RadMom1D_pState_Third_Order W_temp(*this);
  multiplication(W_temp, b);
  return W_temp;
}

inline RadMom1D_pState_Third_Order operator *(const double &b, const RadMom1D_pState_Third_Order &W) {
  RadMom1D_pState_Third_Order W_tmp(W);
  // Use class object U_tmp to access base class member function
  W_tmp.multiplication(b,W_tmp);

  return (W_tmp);
}

//! scalar division
inline RadMom1D_pState_Third_Order RadMom1D_pState_Third_Order :: operator /(const double &b) {
  RadMom1D_pState_Third_Order W_temp(*this);
  division(W_temp, b);
  return W_temp;
}

//! solution state division operator
inline RadMom1D_pState_Third_Order RadMom1D_pState_Third_Order :: operator /(const RadMom1D_pState_Third_Order &W) {
  RadMom1D_pState_Third_Order W_temp(*this);
  division(W_temp, W);
  return W_temp;
}

//! inner product
inline double RadMom1D_pState_Third_Order :: operator *(const RadMom1D_pState_Third_Order &W) {
  return multiplication(*this, W);
}

//! solution state product operator
inline RadMom1D_pState_Third_Order RadMom1D_pState_Third_Order :: operator ^(const RadMom1D_pState_Third_Order &W) {
  RadMom1D_pState_Third_Order W_temp(*this);
  dot_product(W_temp, W);
  return W_temp;
}

/****************************************************************//**
 * RadMom1D_pState_Third_Order -- Assignment operator.
 ********************************************************************/
inline RadMom1D_pState_Third_Order& RadMom1D_pState_Third_Order :: operator =(const RadMom1D_pState_Third_Order &W) {
  equals(*this, W);
  return (*this);
}

/****************************************************************//**
 * RadMom1D_pState_Third_Order -- Shortcut operators.
 ********************************************************************/
inline RadMom1D_pState_Third_Order& RadMom1D_pState_Third_Order :: operator +=(const RadMom1D_pState_Third_Order &W) {
  plus_equal(*this, W);
  return (*this);
}

inline RadMom1D_pState_Third_Order& RadMom1D_pState_Third_Order :: operator -=(const RadMom1D_pState_Third_Order &W) {
  minus_equal(*this, W);
  return (*this);
}

inline RadMom1D_pState_Third_Order& RadMom1D_pState_Third_Order::operator *=(const double &b) {
  times_equal(*this, b);
  return (*this);
}

inline RadMom1D_pState_Third_Order& RadMom1D_pState_Third_Order::operator *=(const RadMom1D_pState_Third_Order &W) {
  times_equal(*this, W);
  return (*this);
}

inline RadMom1D_pState_Third_Order& RadMom1D_pState_Third_Order::operator /=(const double &b) {
  divide_equal(*this, b);
  return (*this);
}

inline RadMom1D_pState_Third_Order& RadMom1D_pState_Third_Order::operator /=(const RadMom1D_pState_Third_Order &W) {
  divide_equal(*this, W);
  return (*this);
}

/********************************************************
 * RadMom1D_pState_Third_Order -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const RadMom1D_pState_Third_Order &W) {
  RadMom1D_pState_Third_Order W_temp;
    W_temp.output(out_file, W);
    return (out_file);
}

inline istream &operator >> (istream &in_file, RadMom1D_pState_Third_Order &W) {
    W.input(in_file, W);
    return (in_file);
}

/****************************************************************//**
 * RadMom1D_cState_Third_Order -- Binary arithmetic operators.
 ********************************************************************/
//! addition
inline RadMom1D_cState_Third_Order RadMom1D_cState_Third_Order :: operator +(const RadMom1D_cState_Third_Order &U) {
  RadMom1D_cState_Third_Order U_temp(*this);
  addition(U_temp, U);
  return U_temp;
}

//! subtraction
inline RadMom1D_cState_Third_Order RadMom1D_cState_Third_Order :: operator -(const RadMom1D_cState_Third_Order &U) {
  RadMom1D_cState_Third_Order U_temp(*this);
  subtraction(U_temp, U);
  return U_temp;
}

//! scalar multiplication
inline RadMom1D_cState_Third_Order RadMom1D_cState_Third_Order :: operator *(const double &b) {
  RadMom1D_cState_Third_Order U_temp(*this);
  multiplication(U_temp, b);
  return U_temp;
}

inline RadMom1D_cState_Third_Order operator *(const double &b, const RadMom1D_cState_Third_Order &U) {
  RadMom1D_cState_Third_Order U_tmp(U);

  // Use class object U_tmp to access base class member function
  U_tmp.multiplication(b,U_tmp);

  return U_tmp;
}

//! scalar division
inline RadMom1D_cState_Third_Order RadMom1D_cState_Third_Order :: operator /(const double &b) {
  RadMom1D_cState_Third_Order U_temp(*this);
  division(U_temp, b);
  return U_temp;
}

//! solution state division operator
inline RadMom1D_cState_Third_Order RadMom1D_cState_Third_Order :: operator /(const RadMom1D_cState_Third_Order &U) {
  RadMom1D_cState_Third_Order U_temp(*this);
  division(U_temp, U);
  return U_temp;
}

//! inner product
inline double RadMom1D_cState_Third_Order :: operator *(const RadMom1D_cState_Third_Order &U) {
  return multiplication(*this, U);
}

//! solution state product operator
inline RadMom1D_cState_Third_Order RadMom1D_cState_Third_Order :: operator ^(const RadMom1D_cState_Third_Order &U) {
  RadMom1D_cState_Third_Order U_temp(*this);
  dot_product(U_temp, U);
  return U_temp;
}

/****************************************************************//**
 * RadMom1D_cState_Third_Order -- Assignment operator.
 ********************************************************************/
inline RadMom1D_cState_Third_Order& RadMom1D_cState_Third_Order :: operator =(const RadMom1D_cState_Third_Order &U) {
  equals(*this, U);
  return (*this);
}

/****************************************************************//**
 * RadMom1D_cState_Third_Order -- Shortcut operators.
 ********************************************************************/
inline RadMom1D_cState_Third_Order& RadMom1D_cState_Third_Order :: operator +=(const RadMom1D_cState_Third_Order &U) {
  plus_equal(*this, U);
  return (*this);
}

inline RadMom1D_cState_Third_Order& RadMom1D_cState_Third_Order :: operator -=(const RadMom1D_cState_Third_Order &U) {
  minus_equal(*this, U);
  return (*this);
}

inline RadMom1D_cState_Third_Order& RadMom1D_cState_Third_Order::operator *=(const double &b) {
  times_equal(*this, b);
  return (*this);
}

inline RadMom1D_cState_Third_Order& RadMom1D_cState_Third_Order::operator *=(const RadMom1D_cState_Third_Order &U) {
  times_equal(*this, U);
  return (*this);
}

inline RadMom1D_cState_Third_Order& RadMom1D_cState_Third_Order::operator /=(const double &b) {
  divide_equal(*this, b);
  return (*this);
}

inline RadMom1D_cState_Third_Order& RadMom1D_cState_Third_Order::operator /=(const RadMom1D_cState_Third_Order &U) {
  divide_equal(*this, U);
  return (*this);
}

/********************************************************
 * RadMom1D_cState_Third_Order -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const RadMom1D_cState_Third_Order &U) {
  RadMom1D_cState_Third_Order U_temp;
    U_temp.output(out_file, U);
    return (out_file);
}

inline istream &operator >> (istream &in_file, RadMom1D_cState_Third_Order &U) {
    U.input(in_file, U);
    return (in_file);
}

/****************************************************************//**
 * Set array pointers to null.
 ********************************************************************/
inline void RadMom1D_pState_Third_Order :: Nullify()
{
#ifndef STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER
  m_values = NULL;
#endif
};

inline void RadMom1D_pState_Third_Order :: Allocate()
{
#ifndef STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER

  // deallocate first
  Deallocate();

  // create the jagged array
  if (NUM_VAR_RADMOM1D_THIRD_ORDER > 0)  { m_values = new double[NUM_VAR_RADMOM1D_THIRD_ORDER]; }

#endif
}


inline void RadMom1D_pState_Third_Order :: Deallocate()
{

#ifndef STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER
  if ( m_values != NULL ) { delete[] m_values; m_values = NULL; }
#endif

}

/****************************************************************//**
 * Set array pointers to null.
 ********************************************************************/
inline void RadMom1D_cState_Third_Order :: Nullify()
{
#ifndef STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER
  m_values = NULL;
#endif
};

inline void RadMom1D_cState_Third_Order :: Allocate()
{
#ifndef STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER

  // deallocate first
  Deallocate();

  // create the jagged array
  if (NUM_VAR_RADMOM1D_THIRD_ORDER > 0)  { m_values = new double[NUM_VAR_RADMOM1D_THIRD_ORDER]; }

#endif
}


inline void RadMom1D_cState_Third_Order :: Deallocate()
{

#ifndef STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER
  if ( m_values != NULL ) { delete[] m_values; m_values = NULL; }
#endif

}

/****************************************************************//**
 * RadMom1D_pState_Third_Order :: Setting up static variables.
 ********************************************************************/
inline void RadMom1D_pState_Third_Order :: SetupStatic( const int &ScatteringFunc,
                                                        const int &i_Moment_Closure,
                                                        const int &i_AbsorptionModel)
{
    // deallocate static vars first, just to be safe
    DeallocateStatic();

    NUM_VAR_RADMOM1D_THIRD_ORDER = STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER;

    RadMom1D_pState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order> :: SetupStatic(ScatteringFunc,
                                                                i_Moment_Closure,
                                                                i_AbsorptionModel);

    if (closure_type != MOMENT_CLOSURE_P3) {
        cout << "Closure type not specified properly !!!!" << endl;
        cout << "Should specify MOMENT_CLOSURE_P3" << endl;
        exit(0);
    }
}

inline void RadMom1D_pState_Third_Order :: DeallocateStatic()
{

}

/****************************************************************//**
 * RadMom1D_cState_Third_Order :: Setting up static variables.
 ********************************************************************/
inline void RadMom1D_cState_Third_Order :: SetupStatic( const int &ScatteringFunc,
                                                        const int &i_Moment_Closure,
                                                        const int &i_AbsorptionModel)
{
    // deallocate static vars first, just to be safe
    DeallocateStatic();

    NUM_VAR_RADMOM1D_THIRD_ORDER = STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER;

    RadMom1D_pState<RadMom1D_cState_Third_Order,
                    RadMom1D_pState_Third_Order> :: SetupStatic(ScatteringFunc,
                                                                i_Moment_Closure,
                                                                i_AbsorptionModel);

    if (closure_type != MOMENT_CLOSURE_P3) {
        cout << "Closure type not specified properly !!!!" << endl;
        cout << "Should specify MOMENT_CLOSURE_P3" << endl;
        exit(0);
    }
}

inline void RadMom1D_cState_Third_Order :: DeallocateStatic()
{

}

/********************************************
 * RadMom1D_pState_Third_Order Value Constructor.        *
 *******************************************/
inline RadMom1D_pState_Third_Order::RadMom1D_pState_Third_Order(const double &Val){
    for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
        m_values[i] = Val;
    }
}

/********************************************
 * RadMom1D_cState_Third_Order Value Constructor.        *
 *******************************************/
inline RadMom1D_cState_Third_Order::RadMom1D_cState_Third_Order(const double &Val){
    for ( int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++ ){
        m_values[i] = Val;
    }
}

/********************************************************
 * RadMom1D_pState_Third_Order::fsca -- Absolute radiative flux.    *
 ********************************************************/
inline double RadMom1D_pState_Third_Order::fsca() const {
  return ( sqrt(f2()) );
}

/********************************************************
 * RadMom1D_pState_Third_Order::f2 -- Square of radiative flux.     *
 ********************************************************/
inline double RadMom1D_pState_Third_Order::f2() const {
  return (sqr(N1x()));
}

/********************************************************
 * RadMom1D_cState_Third_Order::fsca -- Absolute radiative flux.    *
 ********************************************************/
inline double RadMom1D_cState_Third_Order::fsca() const {
  return ( sqrt(f2()) );
}

/********************************************************
 * RadMom1D_cState_Third_Order::f2 -- Square of radiative flux.     *
 ********************************************************/
inline double RadMom1D_cState_Third_Order::f2() const {
    return W().f2();
}

/********************************************************
 * RadMom1D_pState_Third_Order::RadMom1D_pState_Third_Order -- Constructor.       *
 ********************************************************/
inline RadMom1D_pState_Third_Order::RadMom1D_pState_Third_Order(const RadMom1D_cState_Third_Order &U) {
    Nullify();
    Allocate();
    
    if (U.I0() == ZERO){
        m_values[0] = ZERO; 
        m_values[1] = ZERO;
        m_values[2] = ZERO;
        m_values[3] = ZERO;
    } else {
        m_values[0] = U.I0(); 
        m_values[1] = U.I1x()/U.I0();
        m_values[2] = U.I2xx()/U.I0();
        m_values[3] = U.I3xxx()/U.I0();
    }
}

/********************************************************
 * RadMom1D_pState_Third_Order::U -- Conserved solution state.       *  (Start adding closure_type here)
 ********************************************************/
inline RadMom1D_cState_Third_Order RadMom1D_pState_Third_Order::U(void) const {
    return U_from_W(*this);
}

inline RadMom1D_cState_Third_Order RadMom1D_pState_Third_Order::U(const RadMom1D_pState_Third_Order &W) {
    return W.U();
}

inline RadMom1D_cState_Third_Order U(const RadMom1D_pState_Third_Order &W) {
    return W.U();
}


/********************************************************
 * RadMom1D_cState_Third_Order::RadMom1D_cState_Third_Order -- Constructor.     *
 ********************************************************/
inline RadMom1D_cState_Third_Order::RadMom1D_cState_Third_Order(const RadMom1D_pState_Third_Order &W) {
    Nullify();
    Allocate();

    m_values[0] = W.I0();
    m_values[1] = W.I0()*W.N1x();
    m_values[2] = W.I0()*W.N2xx();
    m_values[3] = W.I0()*W.N3xxx();
}

/********************************************************
 * RadMom1D_cState_Third_Order::W -- Primitive solution state.       *
 ********************************************************/
inline RadMom1D_pState_Third_Order RadMom1D_cState_Third_Order::W(void) const {
    return W_from_U(*this);
}

inline RadMom1D_pState_Third_Order RadMom1D_cState_Third_Order::W(const RadMom1D_cState_Third_Order &U) {
    return U.W();
}

inline RadMom1D_pState_Third_Order W(const RadMom1D_cState_Third_Order &U) {
    return U.W();
}

/********************************************************
 * RadMom1D_pState_Third_Order::Fx -- Solution flux (x-direction).   *
 ********************************************************/
inline RadMom1D_cState_Third_Order RadMom1D_pState_Third_Order::Fx(void) const {
    RadMom1D_cState_Third_Order flux;
    
    flux[1] = I0()*N1x();
    flux[2] = I0()*N2xx();
    flux[3] = I0()*N3xxx();
    flux[4] = I0()*r_xxxx();
    
    return flux;
}

inline RadMom1D_cState_Third_Order RadMom1D_pState_Third_Order::Fx(const RadMom1D_pState_Third_Order &W) {
    return W.Fx();
}

inline RadMom1D_cState_Third_Order Fx(const RadMom1D_pState_Third_Order &W) {
    return W.Fx();
}

/************************************************************
 * RadMom1D_pState_Third_Order::lambda_x -- Eigenvalue(s) (x-direction). *
 ************************************************************/
inline RadMom1D_pState_Third_Order RadMom1D_pState_Third_Order::lambda_x(void) const {
    RadMom1D_pState_Third_Order lambda;
    double temp_val_1, temp_val_2;

    temp_val_1 = sqrt((THREE/SEVEN) - (TWO*sqrt(THIRTY)/THIRTY_FIVE));
    temp_val_2 = sqrt((THREE/SEVEN) + (TWO*sqrt(THIRTY)/THIRTY_FIVE));

    lambda[1] = temp_val_1;
    lambda[2] = temp_val_2;
    lambda[3] = -temp_val_1;
    lambda[4] = -temp_val_2;

    return lambda;
}

inline RadMom1D_pState_Third_Order RadMom1D_pState_Third_Order::lambda_x(const RadMom1D_pState_Third_Order &W) {
    return W.lambda_x();
}

inline RadMom1D_pState_Third_Order lambda_x(const RadMom1D_pState_Third_Order &W) {
    return W.lambda_x();
}

//******************************************************************************
// RadMom1D_pState_Third_Order::Setup_Eigenstructure_P3.
//******************************************************************************
inline void RadMom1D_pState_Third_Order :: Setup_Eigenstructure_P3(Eigenstructure_P3 &Eig_P3) const {
    RadMom1D_pState_Third_Order lambda;
    lambda = lambda_x();

    // Compute the eigenvalues of the the flux-Jacobian matrix
    for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++) {
        Eig_P3.lambdas[i] = lambda[i+1];
    }

    // Compute the eigenvectors of the the flux-Jacobian matrix
    // Right eigenvectors
    rc(Eig_P3);
    // Left eigenvectors
    lc(Eig_P3);
}

//******************************************************************************
//RadMom1D_pState_Third_Order::rc -- Conserved right eigenvector (x-direction).                 
//******************************************************************************
inline double RadMom1D_pState_Third_Order::rc(const int &index_U, const int &j) const {
    double rc_vec;

    switch(index_U) {
        case 1 :
            switch(j) {
                case 1 :
                    rc_vec = sqrt(15.0 - 2.0*sqrt(30.0))*(23.0*sqrt(35.0) + 20.0*sqrt(42.0))/21.0;
                    break;
                case 2 :
                    rc_vec = sqrt(15.0 + 2.0*sqrt(30.0))*(23.0*sqrt(35.0) - 20.0*sqrt(42.0))/21.0;
                    break;
                case 3 :
                    rc_vec = -sqrt(15.0 - 2.0*sqrt(30.0))*(23.0*sqrt(35.0) + 20.0*sqrt(42.0))/21.0;
                    break;
                case 4 :
                    rc_vec = -sqrt(15.0 + 2.0*sqrt(30.0))*(23.0*sqrt(35.0) - 20.0*sqrt(42.0))/21.0;
                    break;
                default:
                    rc_vec = ZERO;
                    break;
            }
            break;
        case 2 :
            switch(j) {
                case 1 :
                    rc_vec = FIVE + (TWO*sqrt(THIRTY)/THREE);
                    break;
                case 2 :
                    rc_vec = FIVE - (TWO*sqrt(THIRTY)/THREE);
                    break;
                case 3 :
                    rc_vec = FIVE + (TWO*sqrt(THIRTY)/THREE);
                    break;
                case 4 :
                    rc_vec = FIVE - (TWO*sqrt(THIRTY)/THREE);
                    break;
                default:
                    rc_vec = ZERO;
                    break;
            }
            break;
        case 3 :
            switch(j) {
                case 1 :
                    rc_vec = sqrt(15.0 - 2.0*sqrt(30.0))*(3.0*sqrt(35.0) + 2.0*sqrt(42.0))/21.0;
                    break;
                case 2 :
                    rc_vec = sqrt(15.0 + 2.0*sqrt(30.0))*(3.0*sqrt(35.0) - 2.0*sqrt(42.0))/21.0;
                    break;
                case 3 :
                    rc_vec = -sqrt(15.0 - 2.0*sqrt(30.0))*(3.0*sqrt(35.0) + 2.0*sqrt(42.0))/21.0;
                    break;
                case 4 :
                    rc_vec = -sqrt(15.0 + 2.0*sqrt(30.0))*(3.0*sqrt(35.0) - 2.0*sqrt(42.0))/21.0;
                    break;
                default:
                    rc_vec = ZERO;
                    break;
            }
            break;
        case 4 :
            switch(j) {
                case 1 :
                    rc_vec = ONE;
                    break;
                case 2 :
                    rc_vec = ONE;
                    break;
                case 3 :
                    rc_vec = ONE;
                    break;
                case 4 :
                    rc_vec = ONE;
                    break;
                default:
                    rc_vec = ZERO;
                    break;
            }
            break;
        default:
            cout << "Incorrect value for index_U = " << index_U << endl;
            exit(0);
            break;
    }

    return rc_vec;
}

inline void RadMom1D_pState_Third_Order::rc(Eigenstructure_P3 &Eig_P3) const {
    // Find the vector rc such that: A rc = \lam_i rc
    // where A is the flux Jacobian and \lam_i is the an eigenvalue of the matrix A
    // Right eigenvectors are stored columnwise in the matrix of eigenvectors
    
    for (int i = 1; i <= STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++) {
        for (int j = 1; j <= STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; j++) {
            Eig_P3.rc_vec[i-1][j-1] = rc(i, j);
        }
    }
}

inline void rc(const RadMom1D_pState_Third_Order &W, Eigenstructure_P3 &Eig_P3) {
    W.rc(Eig_P3);
}

//******************************************************************************
//RadMom1D_pState_Third_Order::lc -- Conserved left eigenvector (x-direction).                 
//******************************************************************************
inline double RadMom1D_pState_Third_Order::lc(const int &index_U, const int &j) const {
    double lc_vec;

    switch(index_U) {
        case 1 :
            switch(j) {
                case 1 :
                    lc_vec = sqrt(42.0)*sqrt(15.0 - 2.0*sqrt(30.0))/560.0;
                    break;
                case 2 :
                    lc_vec = sqrt(30.0)/80.0;
                    break;
                case 3 :
                    lc_vec = sqrt(15.0 - 2.0*sqrt(30.0))*(4.0*sqrt(35.0) - 5.0*sqrt(42.0))/560.0;
                    break;
                case 4 :
                    lc_vec = (1.0/4.0) - sqrt(30.0)/16.0;
                    break;
                default:
                    lc_vec = ZERO;
                    break;
            }
            break;
        case 2 :
            switch(j) {
                case 1 :
                    lc_vec = -sqrt(42.0)*sqrt(15.0 + 2.0*sqrt(30.0))/560.0;
                    break;
                case 2 :
                    lc_vec = -sqrt(30.0)/80.0;
                    break;
                case 3 :
                    lc_vec = sqrt(15.0 + 2.0*sqrt(30.0))*(4.0*sqrt(35.0) + 5.0*sqrt(42.0))/560.0;
                    break;
                case 4 :
                    lc_vec = (1.0/4.0) + sqrt(30.0)/16.0;
                    break;
                default:
                    lc_vec = ZERO;
                    break;
            }
            break;
        case 3 :
            switch(j) {
                case 1 :
                    lc_vec = -sqrt(42.0)*sqrt(15.0 - 2.0*sqrt(30.0))/560.0;
                    break;
                case 2 :
                    lc_vec = sqrt(30.0)/80.0;
                    break;
                case 3 :
                    lc_vec = -sqrt(15.0 - 2.0*sqrt(30.0))*(4.0*sqrt(35.0) - 5.0*sqrt(42.0))/560.0;
                    break;
                case 4 :
                    lc_vec = (1.0/4.0) - sqrt(30.0)/16.0;
                    break;
                default:
                    lc_vec = ZERO;
                    break;
            }
            break;
        case 4 :
            switch(j) {
                case 1 :
                    lc_vec = sqrt(42.0)*sqrt(15.0 + 2.0*sqrt(30.0))/560.0;
                    break;
                case 2 :
                    lc_vec = -sqrt(30.0)/80.0;
                    break;
                case 3 :
                    lc_vec = -sqrt(15.0 + 2.0*sqrt(30.0))*(4.0*sqrt(35.0) + 5.0*sqrt(42.0))/560.0;
                    break;
                case 4 :
                    lc_vec = (1.0/4.0) + sqrt(30.0)/16.0;
                    break;
                default:
                    lc_vec = ZERO;
                    break;
            }
            break;
        default:
            cout << "Incorrect value for index_U = " << index_U << endl;
            exit(0);
            break;
    }

    return lc_vec;
}

inline void RadMom1D_pState_Third_Order::lc(Eigenstructure_P3 &Eig_P3) const {
    // Find the vector rc such that: A rc = \lam_i rc
    // where A is the flux Jacobian and \lam_i is the an eigenvalue of the matrix A
    // Left eigenvectors are stored row-wise in the matrix of eigenvectors
    
    for (int i = 1; i <= STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++) {
        for (int j = 1; j <= STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; j++) {
            Eig_P3.lc_vec[i-1][j-1] = lc(i, j);
        }
    }
}

inline void lc(const RadMom1D_pState_Third_Order &W, Eigenstructure_P3 &Eig_P3) {
    W.lc(Eig_P3);
}

//******************************************************************************
//RadMom1D_pState_Third_Order::Compute_Correction_Factor_Roe: compute Correction
// factor for Roe approximate flux function
//******************************************************************************
inline void RadMom1D_pState_Third_Order :: Compute_Correction_Factor_Roe(double *Correction_Factor,
                                                                         const RadMom1D_pState_Third_Order &Wl,
                                                                         const RadMom1D_pState_Third_Order &Wr) const {
    RadMom1D_cState_Third_Order dUrl;
    static double A_Roe[STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER][STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER];

    // Compute the jump in the solution between the left and the right states
    dUrl = Wr.U() - Wl.U();

    Compute_Correction_Factor_Approximate_Jacobian_Roe(A_Roe);

    for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++) {
        Correction_Factor[i] = ZERO;
        for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; j++) {
            Correction_Factor[i] += A_Roe[i][j] * dUrl.m_values[j];
        }
    }
}

inline void RadMom1D_pState_Third_Order :: Compute_Correction_Factor_Approximate_Jacobian_Roe(double A_Roe[][STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER]) const {
    Eigenstructure_P3 Eig_P3;
    RadMom1D_pState_Third_Order Wstar;
    double lambda_val, lc_val, rc_val;
    static double A_Roe_Temp[STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER][STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER];

    Setup_Eigenstructure_P3(Eig_P3);

    for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++) {
        for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; j++) {
            rc_val = Eig_P3.rc_vec[i][j];
            lambda_val = Eig_P3.lambdas[j];
            A_Roe_Temp[i][j] = rc_val*fabs(lambda_val);
        }
    }

    for (int i = 0; i < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; i++) {
        for (int j = 0; j < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; j++) {
            A_Roe[i][j] = ZERO;
            for (int k = 0; k < STATIC_NUM_VAR_RADMOM1D_THIRD_ORDER; k++) {
                lc_val = Eig_P3.lc_vec[k][j];
                A_Roe[i][j] += A_Roe_Temp[i][k]*lc_val;
            }
        }
    }
}

/**********************************************************************
 * RadMom1D_pState_Third_Order::S -- Include all source term vectors and    *
 *                             Jacobians.                             *
**********************************************************************/
inline RadMom1D_cState_Third_Order RadMom1D_cState_Third_Order::S(const Medium1D_State &M ) {
    RadMom1D_cState_Third_Order Source;
    double sigma( M.sigma() );
    
    Source[1] = M.kappa()*(FOUR*PI*M.Ib() - I0());
    Source[2] = -M.beta()*I1x();
    Source[3] = (ONE/THREE)*(M.kappa()*FOUR*PI*M.Ib()+sigma*I0())-M.beta()*I2xx();
    Source[4] = -M.beta()*I3xxx();

    return Source;
}

inline double RadMom1D_cState_Third_Order::Sr(const Medium1D_State &M ) {
    return Srad(*this, M);
}

/*****************************************************************
 * First order piecewise linear solution reconstruction
 ******************************************************************/
inline void RadMom1D_pState_Third_Order::Reconstruct( const RadMom1D_pState_Third_Order &Wc, 
                                                      const RadMom1D_pState_Third_Order &phi, 
                                                      const RadMom1D_pState_Third_Order &dWdx, 
                                                      const double &dX) {
    Solution_Reconstruct( *this, Wc, phi, dWdx, dX);
}

/********************************************************
 * RadMom1D_cState_Third_Order::Fx -- Solution flux (x-direction).   *
 ********************************************************/
inline RadMom1D_cState_Third_Order RadMom1D_cState_Third_Order::Fx(void) const {
    return W().Fx();
}

inline RadMom1D_cState_Third_Order RadMom1D_cState_Third_Order::Fx(const RadMom1D_cState_Third_Order &U) {
    return U.W().Fx();
}

inline RadMom1D_cState_Third_Order Fx(const RadMom1D_cState_Third_Order &U) {
    return U.W().Fx();
}

#endif /* _RADMOM1D_STATE_THIRD_ORDER_INCLUDED  */
