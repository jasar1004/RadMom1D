/////////////////////////////////////////////////////////////////////
///
/// \file Medium1DState.cc
/// 
/// \author J.A.R. Sarr (Jojo)
/// 
/// The radiation state class contains the properties of a gray,
/// absorbing, emitting, anisotropically scattering medium.  This file
/// defines the state class.
///
/////////////////////////////////////////////////////////////////////
#include "Medium1DState.h"
int    Medium1D_State :: NUM_VAR_MEDIUM1D = 0;   // total number of variables
int    Medium1D_State :: Absorb_Type    = MEDIUM1D_ABSORB_GRAY; // absorbsion model
bool   Medium1D_State :: Scatter_Iso    = true;  // isotropic or anisotropic scattering

/////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////

/****************************************************************//**
 * Compute Band Blackbody Intensity 
 ********************************************************************/
void Medium1D_State :: getBlackBody(const double &T, double &ib) {
  if (Absorb_Type == MEDIUM1D_ABSORB_GRAY) {
    ib = BlackBody(T);
  } else {
    cerr << "\nRte1D_State.cc::Gray_Wall(): Invalid value for absorbsion type flag.\n";
    exit(-1);
  } // endif
};

void Medium1D_State :: setBlackBody(const double &T) {
  getBlackBody(T, Ib());
};

/****************************************************************//**
 * Compute values and initialize state.
 ********************************************************************/
void Medium1D_State :: SetInitialValues( const double &Temperature,
                                         const double &AbsorptionCoef,
                                         const double &ScatteringCoef) {
  // Use Gray Gas (ie. constant)
  if (Absorb_Type == MEDIUM1D_ABSORB_GRAY) {
    kappa() = AbsorptionCoef;
    Ib   () = BlackBody(Temperature);
  }

  // error
  else{
    cerr << "Medium1D_State::SetInitialValues() - Invalid flag for Absorbsion model\n";
    exit(1);
  }

  // Medium Temperature
  T() = Temperature;
  
  // scattering coefficient always assumed gray
  sigma() = ScatteringCoef;

}

/****************************************************************//**
 * Setup Static variables.
 ********************************************************************/
void Medium1D_State :: SetupStatic( const int &i_Absorb_Type,
                                    const int &i_Scattering_Type) {
  int Nfreq;
  // set the absorption type flag
  Absorb_Type = i_Absorb_Type;
  
  // GRAY
  if (Absorb_Type == MEDIUM1D_ABSORB_GRAY) {
    Nfreq = 1;
  } 

  // ERROR
  else {
    cerr << "Medium1D_State::SetupState - Invalid flag for gas type\n";
    exit(-1);
  }
  
  // set the isotropic scattering flag
  Scatter_Iso = (i_Scattering_Type == RADIATION_SCATTER_ISO);

  // set number of variables
  NUM_VAR_MEDIUM1D = 2*Nfreq + 2;
  
  // check to make sure we are not asking for too much memory
#ifdef MEDIUM1D_STATIC_NUMBER_OF_BANDS
  if( MEDIUM1D_STATIC_NUMBER_OF_BANDS < NUM_VAR_MEDIUM1D ) {
    cerr <<"\n WARNING USING STATIC MEDIUM1D BUILT WITH "
         << MEDIUM1D_STATIC_NUMBER_OF_BANDS
         <<" VARS PREDEFINED, HOWEVER ASKING FOR "
         << NUM_VAR_MEDIUM1D << endl;
    exit(1); 
  }
#endif

}

/****************************************************************//**
 * Compute NEW values for the medium state dependent on relevant
 * paremters.
 ********************************************************************/
void Medium1D_State :: SetState( const double &Temperature ) {
  // Use Gray Gas (ie. constant)
  if (Absorb_Type == MEDIUM1D_ABSORB_GRAY) {
    //kappa() = AbsorptionCoef; // <-- it is a treated as a specified constant
    Ib() = BlackBody(Temperature);
  }

  // error
  else{
    cerr << "Medium1D_State::SetInitialValues() - Invalid flag for Absorbsion model\n";
    exit(1);
  }

  // Medium Temperature
  T() = Temperature;
  
  // Scattering coefficient 

  // scattering coefficient always assumed gray
  // sigma() = ScatteringCoef; // <-- it is a treated specified constant
  
}

/********************************************************
 * Blackbody intensity                                  *
 ********************************************************/
double Medium1D_State :: BlackBody(const double T) {
  return STEFFAN_BOLTZMANN*pow(T,4.0)/PI;
}
