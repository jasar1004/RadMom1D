/*!\file RadMom1DState.h
  \brief Header file defining 1D RadMom Solution State Classes. */

#ifndef _RADMOM1D_STATE_INCLUDED
#define _RADMOM1D_STATE_INCLUDED

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

#ifndef _CFD_INCLUDED
#include "./CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MEDIUM_1D_STATE_INCLUDED
#include "./Medium1DState/Medium1DState.h"
#endif //_MEDIUM_1D_STATE_INCLUDED

// Enviroment flag for cfrte1D root directory path
#define PATHVAR "cfrte1D_Path"

/* Define the classes. */
template<class cState, class pState>
class RadMom1D_cState;

template<class cState, class pState>
class RadMom1D_pState;

/*!
 * Class: RadMom1D_pState
 *
 * @brief Base class for primitive variable solution state class .
 *
 * \verbatim
 * Member functions
 *     k        -- Return Boltzmann constant.
 *     h        -- Return Planck's constant.
 *     c        -- Return speed of light.
 *     sig      -- Return addition of absorption and scattering coefficients.
 *     a        -- Return 8*pi^5*k^4/(15*h^3*c^3).
 *     Sr        -- Return source term.
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

template<class cState, class pState>
class RadMom1D_pState {
protected:

private:
public:
  static int              closure_type; //!< Identifier for either M1 or P1 moment closure.
  static int              Absorption_Model;
  static int              Scattering_Func;
  static double                    c; //!< Speed of light.
  static double                    a; //!< Blackbody source constant.

  // Planck constants
  static double                   C1; // [W/(m2 ster m-4)]

  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  RadMom1D_pState(void) {}

  /* Destructor. */
  ~RadMom1D_pState(void) { }
  // Use automatically generated destructor.
  //@}

  //! Copy operator.
  pState Copy(const pState &W) {
      pState W_temp;
      for ( int i = 0; i < W.NumVar(); i++ ){
           W_temp.m_values[i] = W.m_values[i];
      }
      return W_temp;
  }

  void Copy_to_W(pState &W_copy, const pState &W) {
      for ( int i = 0; i < W.NumVar(); i++ ){
           W_copy.m_values[i] = W.m_values[i];
      }
  }

  //! Set array pointers to null
  void Nullify(pState &W);

  //! Vacuum operator.
  void Vacuum(pState &W) {
      for ( int i = 0; i < W.NumVar(); i++ ){
           W.m_values[i] = ZERO;
      }
  }
  //! Ones operator.
  void Ones(pState &W) {
      for ( int i = 0; i < W.NumVar(); i++ ){
           W.m_values[i] = ONE;
      }
  }

  //@}

  //@{ @name Conserved solution state.
  double U_from_W(const pState &W, const int &index_U) const;
  cState U_from_W(const pState &W) const;

  // friend cState U(const pState &W);
  //@}

  //@{ @name First order linear solution reconstruction
  void Solution_Reconstruct( pState &Wrecon,
                             const pState &Wc,
                             const pState &phi,
                             const pState &dWdx,
                             const double &dX);

    //! @name Allocators and deallocators
  //@{

  //! Setup state static variables
  //! memory allocation / deallocation for the static arrays
  static void SetupStatic( const int &ScatteringFunc,
                           const int &i_Moment_Closure,
                           const int &i_AbsorptionModel);
  static void DeallocateStatic();

  //@{ @name Binary arithmetic operators.
  void addition(pState &W_this, const pState &W2);
  void subtraction(pState &W_this, const pState &W2);
  void multiplication(pState &W_this, const double &b);
  void multiplication(const double &b, pState &W_this);
  void division(pState &W_this, const double &b);
  void division(pState &W_this, const pState &W_denom);
  double multiplication(pState &W_this, const pState &W2);
  void dot_product(pState &W_this, const pState &W2);

  /* Assignment operator. */
  void equals(pState &W_this, const pState &W);

  //@{ @name Shortcut arithmetic operators.
  void plus_equal(pState &W_this, const pState &W2);
  void minus_equal(pState &W_this, const pState &W2);
  void times_equal(pState &W_this, const double &b);
  void times_equal(pState &W_this, const pState &W2);
  void divide_equal(pState &W_this, const double &b);
  void divide_equal(pState &W_this, const pState &W2);
  //@}

  //@{ @name Input-output operators.
  void output(ostream &out_file, const pState &W);
  void input(istream &in_file,  pState &W);
  //@}

  //! Index operator.
  double& W_index(pState &W, int index)       {
      assert( index >= 1 && index <= W.NumVar() );
      return W.m_values[index-1];
  }
  const double& W_index(const pState &W, int index) const {
      assert( index >= 1 && index <= W.NumVar() );
      return W.m_values[index-1];
  }

  // Routines for computing approximate Roe matrix
  pState U_to_W_Roe(const cState &U) const;

  cState W_Roe_to_U(const pState &W_Roe) const;

};

 // end of RadMom1D_pState class

/*!
 * Class: RadMom1D_cState
 *
 * @brief Base class for conserved variable solution state class .
 *
 * \verbatim
 * Member functions
 *     k        -- Return Boltzmann constant.
 *     h        -- Return Planck's constant.
 *     c        -- Return speed of light.
 *     sig      -- Return addition of absorption and scattering coefficients.
 *     a        -- Return 8*pi^5*k^4/(15*h^3*c^3).
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

template<class cState, class pState>
class RadMom1D_cState {
protected:
private:
public:
  static int              closure_type; //!< Identifier for either M1 or P1 moment closure.
  static int          Absorption_Model;
  static int          Scattering_Func;
  static double                      c; //!< Speed of light.
  static double                      a; //!< Blackbody source constant.

  // Planck constants
  static double                     C1; // [W/(m2 ster m-4)]static

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  RadMom1D_cState(void) { }

  /* Destructor. */
  ~RadMom1D_cState(void) { }
  // Use automatically generated destructor.
  //@}

  //! Copy operator.
  cState Copy(const cState &U) {
      cState U_temp;
      for ( int i = 0; i < U.NumVar(); i++ ){
          U_temp.m_values[i] = U.m_values[i];
      }
      return U_temp;
  }

  void Copy_to_U(cState &U_copied, const cState &U) const {
      for ( int i = 0; i < U.NumVar(); i++ ){
           U_copied.m_values[i] = U.m_values[i];
      }
  }

  //! Set array pointers to null
  void Nullify(pState &W);

  //! Vacuum operator.
  void Vacuum(cState &U) {
      for ( int i = 0; i < U.NumVar(); i++ ){
           U.m_values[i] = ZERO;
      }
  }
  //! Ones operator.
  void Ones(cState &U) {
      for ( int i = 0; i < U.NumVar(); i++ ){
           U.m_values[i] = ONE;
      }
  }

    //@{ @name Primitive solution state.
  pState W_from_U(const cState &U) const ;
  // friend pState W <cState, pState>(const cState &U);
  //@}

  //@{ @name Include all source vectors and Jacobians.
  double Srad(const cState &U, const Medium1D_State &M) const;

  //! Setup state static variables
  //! memory allocation / deallocation for the static arrays
  static void SetupStatic( const int &ScatteringFunc,
                           const int &i_Moment_Closure,
                           const int &i_AbsorptionModel);
  static void DeallocateStatic();

  //@{ @name Binary arithmetic operators.
  void addition(cState &U_this, const cState &U2);
  void subtraction(cState &U_this, const cState &U2);
  void multiplication(cState &U_this, const double &b);
  void multiplication(const double &b, cState &U_this);
  void division(cState &U_this, const double &b);
  void division(cState &U_this, const cState &U_denom);
  double multiplication(cState &U_this, const cState &U2);
  void dot_product(cState &U_this, const cState &U2);

  /* Assignment operator. */
  void equals(cState &U_this, const cState &U);

  //@{ @name Shortcut arithmetic operators.
  void plus_equal(cState &U_this, const cState &U2);
  void minus_equal(cState &U_this, const cState &U2);
  void times_equal(cState &U_this, const double &b);
  void times_equal(cState &U_this, const cState &U2);
  void divide_equal(cState &U_this, const double &b);
  void divide_equal(cState &U_this, const cState &U2);
  //@}

  //@{ @name Input-output operators.
  void output(ostream &out_file, const cState &U);
  void input(istream &in_file,  cState &U);
  //@}

  //! Index operator.
  double& U_index(cState &U, int index)       {
      assert( index >= 1 && index <= U.NumVar() );
      return U.m_values[index-1];
  }
  const double& U_index(const cState &U, int index) const {
      assert( index >= 1 && index <= U.NumVar() );
      return U.m_values[index-1];
  }
};

// End of RadMom1D_cState class

/****************************************************************//**
 * RadMom1D_pState -- Binary arithmetic operators.
 ********************************************************************/
//! addition
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: addition(pState &W_this, const pState &W2) {
  plus_equal(W_this, W2);
}

//! subtraction
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: subtraction(pState &W_this, const pState &W2) {
  minus_equal(W_this, W2);
}

//! scalar multiplication
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: multiplication(pState &W_this, const double &b) {
    times_equal(W_this, b);
}

template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: multiplication(const double &b, pState &W_this) {
    times_equal(W_this, b);
}

//! scalar division
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: division(pState &W_this, const double &b) {
  divide_equal(W_this, b);
}

//! solution state division operator
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: division(pState &W_this, const pState &W_denom) {
    divide_equal(W_this, W_denom);
}

//! inner product
template<class cState, class pState>
inline double RadMom1D_pState<cState, pState> :: multiplication(pState &W1, const pState &W2) {
  double sum=0.0;

  for ( int i = 0; i < W1.NumVar(); i++ ) {
    sum += W1.m_values[i]*W2.m_values[i];
  }
  return sum;
}

//! solution state product operator
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: dot_product(pState &W_this, const pState &W2) {
  times_equal(W_this, W2);
}

/****************************************************************//**
 * RadMom1D_pState -- Assignment operator.
 ********************************************************************/
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: equals(pState &W_this, const pState &W) {
    Copy_to_W(W_this, W);
}

/****************************************************************//**
 * RadMom1D_pState -- Shortcut operators.
 ********************************************************************/
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: plus_equal(pState &W_this, const pState &W2) {
    for ( int i = 0; i < W_this.NumVar(); i++ ) {
        W_this.m_values[i] += W2.m_values[i];
    }
}

template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: minus_equal(pState &W_this, const pState &W2) {
    for ( int i = 0; i < W_this.NumVar(); i++ ) {
        W_this.m_values[i] -= W2.m_values[i];
    }
}

template<class cState, class pState>
inline void RadMom1D_pState<cState, pState>::times_equal(pState &W_this, const double &b) {
    for ( int i = 0; i < W_this.NumVar(); i++ ) {
        W_this.m_values[i] *= b;
    }
}

template<class cState, class pState>
inline void RadMom1D_pState<cState, pState>::times_equal(pState &W_this, const pState &W2) {
    for ( int i = 0; i < W_this.NumVar(); i++ ) {
        W_this.m_values[i] *= W2.m_values[i];
    }
}

template<class cState, class pState>
inline void RadMom1D_pState<cState, pState>::divide_equal(pState &W_this, const double &b) {
    for ( int i = 0; i < W_this.NumVar(); i++ ) {
        W_this.m_values[i] /= b;
    }
}

template<class cState, class pState>
inline void RadMom1D_pState<cState, pState>::divide_equal(pState &W_this, const pState &W2) {
    for ( int i = 0; i < W_this.NumVar(); i++ ) {
        W_this.m_values[i] /= W2.m_values[i];
    }
}

/****************************************************************//**
 * RadMom1D_cState -- Binary arithmetic operators.
 ********************************************************************/
//! addition
template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: addition(cState &U_this, const cState &U2) {
  plus_equal(U_this, U2);
}

//! subtraction
template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: subtraction(cState &U_this, const cState &U2) {
  minus_equal(U_this, U2);
}

//! scalar multiplication
template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: multiplication(cState &U_this, const double &b) {
    times_equal(U_this, b);
}

template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: multiplication(const double &b, cState &U_this) {
    times_equal(U_this, b);
}

//! scalar division
template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: division(cState &U_this, const double &b) {
  divide_equal(U_this, b);
}

//! solution state division operator
template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: division(cState &U_this, const cState &U_denom) {
    divide_equal(U_this, U_denom);
}

//! inner product
template<class cState, class pState>
inline double RadMom1D_cState<cState, pState> :: multiplication(cState &U1, const cState &U2) {
  double sum=0.0;

  for ( int i = 0; i < U1.NumVar(); i++ ) {
    sum += U1.m_values[i]*U2.m_values[i];
  }
  return sum;
}

//! solution state product operator
template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: dot_product(cState &U_this, const cState &U2) {
  times_equal(U_this, U2);
}

/****************************************************************//**
 * RadMom1D_cState -- Assignment operator.
 ********************************************************************/
template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: equals(cState &U_this, const cState &U) {
    Copy_to_U(U_this, U);
}

/****************************************************************//**
 * RadMom1D_cState -- Shortcut operators.
 ********************************************************************/
template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: plus_equal(cState &U_this, const cState &U2) {
    for ( int i = 0; i < U_this.NumVar(); i++ ) {
        U_this.m_values[i] += U2.m_values[i];
    }
}

template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: minus_equal(cState &U_this, const cState &U2) {
    for ( int i = 0; i < U_this.NumVar(); i++ ) {
        U_this.m_values[i] -= U2.m_values[i];
    }
}

template<class cState, class pState>
inline void RadMom1D_cState<cState, pState>::times_equal(cState &U_this, const double &b) {
    for ( int i = 0; i < U_this.NumVar(); i++ ) {
        U_this.m_values[i] *= b;
    }
}

template<class cState, class pState>
inline void RadMom1D_cState<cState, pState>::times_equal(cState &U_this, const cState &U2) {
    for ( int i = 0; i < U_this.NumVar(); i++ ) {
        U_this.m_values[i] *= U2.m_values[i];
    }
}

template<class cState, class pState>
inline void RadMom1D_cState<cState, pState>::divide_equal(cState &U_this, const double &b) {
    for ( int i = 0; i < U_this.NumVar(); i++ ) {
        U_this.m_values[i] /= b;
    }
}

template<class cState, class pState>
inline void RadMom1D_cState<cState, pState>::divide_equal(cState &U_this, const cState &U2) {
    for ( int i = 0; i < U_this.NumVar(); i++ ) {
        U_this.m_values[i] /= U2.m_values[i];
    }
}

/********************************************************
 * RadMom1D_pState -- Input-output operators.            *
 ********************************************************/
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: output(ostream &out_file, const pState &W) {
    out_file.setf(ios::scientific);

    out_file << " ";
    for (int i = 0; i < W.NumVar(); i++) {
        out_file << W[i+1];
        if (i < W.NumVar() - 1) {
            out_file << " ";
        }
    }
    out_file << "\n";

    out_file.unsetf(ios::scientific);
}

template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: input(istream &in_file, pState &W) {
    in_file.setf(ios::skipws);

    for (int i = 0; i < W.NumVar(); i++) {
        in_file >> W[i+1];
    }

    in_file.unsetf(ios::skipws);
}

/********************************************************
 * RadMom1D_cState -- Input-output operators.            *
 ********************************************************/
template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: output(ostream &out_file, const cState &U) {
    out_file.setf(ios::scientific);

    out_file << " ";
    for (int i = 0; i < U.NumVar(); i++) {
        out_file << U[i+1];
        if (i < U.NumVar() - 1) {
            out_file << " ";
        }
    }
    out_file << "\n";

    out_file.unsetf(ios::scientific);
}

template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: input(istream &in_file, cState &U) {
    in_file.setf(ios::skipws);

    for (int i = 0; i < U.NumVar(); i++) {
        in_file >> U[i+1];
    }

    in_file.unsetf(ios::skipws);
}

/****************************************************************//**
 * RadMom1D_pState :: Setting up static variables.
 ********************************************************************/
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: SetupStatic( const int &ScatteringFunc,
                                                            const int &i_Moment_Closure,
                                                            const int &i_AbsorptionModel)
{
    closure_type = i_Moment_Closure;
    Absorption_Model = i_AbsorptionModel;
    Scattering_Func = ScatteringFunc;

    if (closure_type != MOMENT_CLOSURE_P1 &&
        closure_type != MOMENT_CLOSURE_M1 &&
        closure_type != MOMENT_CLOSURE_P3) {
        cout << "Closure type not specified properly !!!!" << endl;
        cout << "Should specify MOMENT_CLOSURE_M1 or MOMENT_CLOSURE_P1 or MOMENT_CLOSURE_P3" << endl;
        exit(0);
    }
}

template<class cState, class pState>
inline void RadMom1D_pState<cState, pState> :: DeallocateStatic()
{

}

/****************************************************************//**
 * RadMom1D_cState :: Setting up static variables.
 ********************************************************************/
template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: SetupStatic( const int &ScatteringFunc,
                                                        const int &i_Moment_Closure,
                                                        const int &i_AbsorptionModel)
{
    // deallocate static vars first, just to be safe
    DeallocateStatic();

    closure_type = i_Moment_Closure;
    Absorption_Model = i_AbsorptionModel;
    Scattering_Func = ScatteringFunc;

    if (closure_type != MOMENT_CLOSURE_P1 &&
        closure_type != MOMENT_CLOSURE_M1 &&
        closure_type != MOMENT_CLOSURE_P3) {
        cout << "Closure type not specified properly !!!!" << endl;
        cout << "Should specify MOMENT_CLOSURE_M1 or MOMENT_CLOSURE_P1 or MOMENT_CLOSURE_P3" << endl;
        exit(0);
    }
}

template<class cState, class pState>
inline void RadMom1D_cState<cState, pState> :: DeallocateStatic()
{

}

/********************************************************
 * RadMom1D_pState::U -- Conserved solution state.       *  (Start adding closure_type here)
 ********************************************************/
template<class cState, class pState>
inline double RadMom1D_pState<cState, pState> :: U_from_W(const pState &W, const int &index_U) const {
    double U_val;
    switch (index_U) {
        case 1:
            U_val = W.m_values[0];
            break;
        case 2:
            U_val = W.m_values[0]*W.m_values[index_U-1];
            break;
        default:
            cout << "Incorrect value for index_U = " << index_U << endl;
            exit(0);
            break;
    };

    return U_val;
}

template<class cState, class pState>
inline cState RadMom1D_pState<cState, pState>::U_from_W(const pState &W) const {
    cState U_temp;

    U_temp[1] = W.m_values[0];
    for (int  i = 1; i < W.NumVar(); i++) {
        U_temp[i+1] = W.m_values[0]*W.m_values[i];
    }

    return U_temp;
}

template<class cState, class pState>
inline cState U(const pState &W) {
    return W.U();
}

/********************************************************
 * RadMom1D_cState::W -- Primitive solution state.       *
 ********************************************************/
template<class cState, class pState>
inline pState RadMom1D_cState<cState, pState>::W_from_U(const cState &U) const {
    pState W_temp;
    W_temp[1] = U.m_values[0];
    for (int  i = 1; i < U.NumVar(); i++) {
        W_temp[i+1] = U.m_values[i]/U.m_values[0];
    }
    return W_temp;
}

template<class cState, class pState>
inline pState W(const cState &U) {
    return U.W();
}

/**********************************************************************
 * RadMom1D_pState::S -- Include all source term vectors and    *
 *                             Jacobians.                             *
 * Regular Source Term
**********************************************************************/
template<class cState, class pState>
inline double RadMom1D_cState<cState, pState>::Srad(const cState &U, const Medium1D_State &M ) const {
    double Source_term;
    Source_term = M.kappa()*(FOUR*PI*M.Ib() - U.I0());

    return -Source_term;
}

/*****************************************************************
 * First order piecewise linear solution reconstruction
 ******************************************************************/
template<class cState, class pState>
inline void RadMom1D_pState<cState, pState>::Solution_Reconstruct( pState &Wrecon,
                                                          const pState &Wc,
                                                          const pState &phi,
                                                          const pState &dWdx,
                                                          const double &dX) {
    for (int  i = 0; i < Wc.NumVar(); i++) {
        Wrecon.m_values[i] = Wc.m_values[i] + phi.m_values[i]*dWdx.m_values[i]*dX;
    }
}
#endif /* _RADMOM1D_STATE_INCLUDED  */
