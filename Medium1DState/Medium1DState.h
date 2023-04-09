/////////////////////////////////////////////////////////////////////
///
/// \file Medium1DState.h
/// 
/// \author J.A.R. Sarr (Jojo)
/// 
/// The medium state class contains the properties  
/// of a gray, absorbing, emitting, anisotropically 
/// scattering medium.   This file defines the state
/// class.                                          
/////////////////////////////////////////////////////////////////////
#ifndef _MEDIUM1D_STATE_INCLUDED
#define _MEDIUM1D_STATE_INCLUDED 

/////////////////////////////////////////////////////////////////////
// INCLUDES
/////////////////////////////////////////////////////////////////////
// Required C++ libraries
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>
using namespace std;

// #include "./RadiatingGas.h"
#include "../Math.h"

#define SPEED_OF_LIGHT      3.00e08
#define RADIATION_CONSTANT  7.56576911e-16
#define BLACKBODY_CONSTANT  5.67037321e-8
#define PLANCK_CONSTANT     3.7418E-16
#define EPSILON_RADMOM      1.00e-5
#define BOLTZMANN           1.380658e-23
#define STEFFAN_BOLTZMANN   5.670e-08

#define TOLER_RADIATIVE_DENSITY   1e-08

/////////////////////////////////////////////////////////////////////
// DEFINES
/////////////////////////////////////////////////////////////////////
//! Absorbsion model type
enum Gas_Models { MEDIUM1D_ABSORB_GRAY,
                  MEDIUM1D_ABSORB_SNBCK
                };

//! Scatter Phase Functions
enum Scatter_Models { RADIATION_SCATTER_ISO,  // isotropic scattering
                      RADIATION_SCATTER_LINEAR,  // linear anisotropic scattering
                    };


//! Moment closure model type
enum Closure_Types { MOMENT_CLOSURE_M1,
                     MOMENT_CLOSURE_P1,
                     MOMENT_CLOSURE_M2,
                     MOMENT_CLOSURE_M2_PROJECTION,
                     MOMENT_CLOSURE_P3
                   };

//! Set fixed static number parameters
#define MEDIUM1D_STATIC_NUMBER_OF_BANDS  4 

/////////////////////////////////////////////////////////////////////
// CLASS DEFINITIONS
/////////////////////////////////////////////////////////////////////

/**
 * \class Medium_State
 *
 * Gas state class definition for a absorbing, scattering 
 * participating medium.
 *
 */
class Medium1D_State {

 public:

  //!< state array
#ifdef MEDIUM1D_STATIC_NUMBER_OF_BANDS
  double x[MEDIUM1D_STATIC_NUMBER_OF_BANDS];
#else 
  double* x;
#endif

  //! @name Static objects
  //@{ 
  static int NUM_VAR_MEDIUM1D; //!< the total number of variables in the state
  static int Absorb_Type;      //!< flag for absorption model
  static bool Scatter_Iso;     //!< true->isotropic scattering, false->anisotropic
  //@}

  //! @name accessors for state array
  //@{ 
  //! blackbody intentsity [W/m^2 or W/(m^2 cm)]
  double        Ib() const { return  x[0]; };
  double&       Ib()       { return  x[0]; };

  //! absorbsion coefficient [m^-1]
  double        kappa() const { return  x[1]; };
  double&       kappa()       { return  x[1]; };

  //! scattering coefficient [m^-1]
  double  sigma(void) const { return x[2]; };
  double& sigma(void)       { return x[2]; };
  //@}

   
  //! Medium Temperature [K]
  double  T(void) const { return x[3]; };
  double&  T(void)      { return x[3]; };
  
  //-----------------------------------------------------------------

  //! @name Constructors/Destructors
  //@{

  Medium1D_State() { Nullify(); Allocate(); Zero(); }

  Medium1D_State( const Medium1D_State &U )
  { 
    Nullify();
    Allocate(); 
    if( this != &U) Copy(U); 
  }

  ~Medium1D_State() { Deallocate(); }

  //! Copy function
  void Copy( const Medium1D_State &U )
  { for ( int i=0; i<NUM_VAR_MEDIUM1D; i++ ) x[i] = U.x[i];}

  //! Set array pointers to null
  void Nullify(void);

  //! Zero operator.
  void Zero() { for(int i=0; i<NUM_VAR_MEDIUM1D; i++) x[i] = 0.0; }

  //! Initialer
  void SetInitialValues( const double &Temperature,
                         const double &AbsorptionCoef,
                         const double &ScatteringCoef);

  //! Return the number of variables.
  static double NumVar( void ) { return NUM_VAR_MEDIUM1D; };

  //@}

  //-----------------------------------------------------------------
  
  //! @name Allocators and deallocators
  //@{ 
  void Allocate();
  void Deallocate();
  //@}

  //! @name Static Functions
  //@{ 

  //! Setup function
  static void SetupStatic( const int &i_Absorb_Type, 
                           const int &i_Scattering_Type);
  //@}

  //-----------------------------------------------------------------

  //! @name State functions.
  //@{ 

  //! Return extinction coefficient
  double beta() const { return (kappa() + sigma()); }

  //! Compute Band Blackbody Intensity
  void getBlackBody(const double &T, double &ib);
  void setBlackBody(const double &T);
  double BlackBody(const double T );

  //! Compute mean absorbsion coefficients
  double PlanckMean(void) const;

  //! Compute a new state dependant upon gas state
  void SetState( const double &Temperature );

  //! print variable list
  static void outputTecplotVarList(ostream &out,
				   const string &prefix = "");
  //@}

  //! deallocate all static variables
  static void DeallocateStatic() {};

  //-----------------------------------------------------------------

  //! @name Operators Overloading
  //@{ 

  //! Assignment operator.
  Medium1D_State& operator =(const Medium1D_State &U) {
    if( this != &U) Copy(U);
    return (*this);
  }

  //! Index operator.
        double& operator[](int index)       { return x[index-1]; };
  const double& operator[](int index) const { return x[index-1]; };

  //! Binary arimatic operators
  Medium1D_State operator -(const Medium1D_State &U) const;
  Medium1D_State operator +(const Medium1D_State &U) const;
  Medium1D_State operator *(const double &a) const;
  friend Medium1D_State operator *(const double &a, const Medium1D_State &U);
  Medium1D_State operator /(const double &a) const;

 //! Shortcut operators
  Medium1D_State& operator -=(const Medium1D_State &U);
  Medium1D_State& operator +=(const Medium1D_State &U);
  Medium1D_State& operator *=(const double &a);
  Medium1D_State& operator /=(const double &a);

  //! Input-output operators.
  void outputTecplot(ostream &out);
  friend ostream& operator << (ostream &out_file, const Medium1D_State &U);
  friend istream& operator >> (istream &in_file,  Medium1D_State &U);

 //@}
};

/////////////////////////////////////////////////////////////////////
// FUNCTION DEFINITIONS
/////////////////////////////////////////////////////////////////////

/****************************************************************//**
 * Set array pointers to null.
 ********************************************************************/
inline void Medium1D_State :: Nullify() 
{ 
#ifndef MEDIUM1D_STATIC_NUMBER_OF_BANDS
  x = NULL;
#endif
};


/****************************************************************//**
 * Array allocator and deallocator for dynamic arrays.
 ********************************************************************/
inline void Medium1D_State :: Allocate()
{
#ifndef MEDIUM1D_STATIC_NUMBER_OF_BANDS

  // deallocate first
  Deallocate();

  // create the arrays
  if (NUM_VAR_MEDIUM1D>0) x = new double[NUM_VAR_MEDIUM1D];
#endif
}

inline void Medium1D_State :: Deallocate() {
#ifndef MEDIUM1D_STATIC_NUMBER_OF_BANDS
  if (x != NULL) { delete[] x; x = NULL; }
#endif
}

/****************************************************************//**
 * Binary arithmetic operators
 ********************************************************************/
inline Medium1D_State Medium1D_State::operator -(const Medium1D_State &U) const{
  Medium1D_State Temp(*this);
  Temp -= U;
  return Temp;
}

inline Medium1D_State Medium1D_State::operator +(const Medium1D_State &U) const{
  Medium1D_State Temp(*this);
  Temp += U;
  return Temp;
}

inline Medium1D_State Medium1D_State::operator *(const double &a) const{
  Medium1D_State Temp(*this);
  Temp *= a;
  return Temp;
}

inline Medium1D_State operator *(const double &a, const Medium1D_State &U) {
  Medium1D_State Temp(U);
  Temp *= a;
  return Temp;
}

inline Medium1D_State Medium1D_State::operator /(const double &a) const{
  Medium1D_State Temp(*this);
  Temp /= a;
  return Temp;
}


/****************************************************************//**
 * Shortcut arithmetic operators
 ********************************************************************/
inline Medium1D_State& Medium1D_State::operator +=(const Medium1D_State &U){
  for(int i=0; i<NUM_VAR_MEDIUM1D; i++) { x[i] += U.x[i]; }
  return (*this);
}

inline Medium1D_State& Medium1D_State::operator -=(const Medium1D_State &U) {
  for(int i=0; i<NUM_VAR_MEDIUM1D; i++) { x[i] -= U.x[i]; }
  return (*this);
}

inline Medium1D_State& Medium1D_State::operator *=(const double &a) {
  for(int i=0; i<NUM_VAR_MEDIUM1D; i++) x[i] *= a;  
  return (*this);
}

inline Medium1D_State& Medium1D_State::operator /=(const double &a) {
  for(int i=0; i<NUM_VAR_MEDIUM1D; i++) x[i] /= a;  
  return (*this);
}

/****************************************************************//**
 * Input/Output operators
 ********************************************************************/
inline ostream& operator << (ostream &out_file, const Medium1D_State &U)
{
  //out_file.precision(10);
  out_file.setf(ios::scientific);
  for( int i=0; i<U.NUM_VAR_MEDIUM1D; i++) out_file<<" "<<U.x[i];
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream& operator >> (istream &in_file,  Medium1D_State &U)
{
  in_file.setf(ios::skipws);
  for( int i=0; i<U.NUM_VAR_MEDIUM1D; i++) in_file >> U.x[i];
  in_file.unsetf(ios::skipws);
  return (in_file);
}

#endif // _MEDIUM1D_STATE_INCLUDED 
