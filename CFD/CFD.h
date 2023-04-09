/*!file CFD.h
  \brief Header file defining CFD subroutines and macros. */

#ifndef _CFD_INCLUDED
#define _CFD_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include <ctime>

using namespace std;

#ifndef _MATH_MACROS_INCLUDED
#include "../Math.h"
#endif // _MATH_MACROS_INCLUDED

/**********************************************************************
 * CFD -- General Purpose Command Codes.                              *
 **********************************************************************/

#define	OFF                           0
#define	ON                            1

/**********************************************************************
 * CFD -- Input parameter command codes.                              *
 **********************************************************************/

#define EXECUTE_CODE                                     10000
#define TERMINATE_CODE                                   10001
#define CONTINUE_CODE                                    10002
#define COMMENT_CODE                                     10003
#define	INVALID_INPUT_CODE                               10004
#define	INVALID_INPUT_VALUE                              10005
#define RUN_TESTS_CODE                                   10006

#define WRITE_OUTPUT_CODE                                10007
#define WRITE_OUTPUT_CELLS_CODE                          10008
#define WRITE_OUTPUT_NODES_CODE                          10009
#define WRITE_RESTART_CODE                               10010

/**********************************************************************
 * CFD -- I/O Types.                                                  *
 **********************************************************************/

#define	IO_TECPLOT       1

/**********************************************************************
 * CFD -- Grid Types.                                                 *
 **********************************************************************/

/**********************************************************************
 * CFD -- Boundary Condition Types.                                   *
 **********************************************************************/

//----DEFAULT NAMES----//                          Core | Scope

//-- Radiation Boundary Conditions
#define BC_NONE                         10
#define BC_FIXED                        11
#define BC_CONSTANT_EXTRAPOLATION       12
#define BC_LINEAR_EXTRAPOLATION         13
#define BC_CHARACTERISTIC               14
#define BC_REFLECTION                   15
#define BC_GRAY_WALL                    16
#define BC_MARSHAK                      17
#define BC_PARTIAL_FLUX                 18
#define BC_PARTIAL_MOMENTS              19
#define BC_INLET_OUTLET                 20
#define BC_AXISYMMETRIC                 21
#define BC_RADIATING_BEAM               22

/**********************************************************************
 * CFD - BC OPTIONS                                                   *
 **********************************************************************/
// Type of boundary conditions enforcement
#define WEAK_BCS                          100
#define STRONG_BCS                        110

/**********************************************************************
 * CFD -- Initial Condition Types.                                    *
 **********************************************************************/

#define IC_NONE                        -3
#define IC_RESTART                     -2
#define IC_CONSTANT                    -1
#define	IC_UNIFORM                      0
#define IC_PARALLEL_PLATES              1
#define RADMOM_ICS_INTENSITY            2
#define RADMOM_ICS_TEMPERATURE          3

/********************************************************
 * CFD -- Time Integration (Time-Stepping) Types.       *
 ********************************************************/
#define	TIME_STEPPING_EXPLICIT_EULER                                 0

/**********************************************************************
 * CFD -- Reconstruction Types.                                       *
 **********************************************************************/

#define	RECONSTRUCTION_GREEN_GAUSS                             1
#define	RECONSTRUCTION_LEAST_SQUARES                           2

/**********************************************************************
 * CFD -- Limiter Types.                                              *
 **********************************************************************/

#define	LIMITER_ONE                                   -1
#define	LIMITER_UNLIMITED                             -1
#define	LIMITER_ZERO                                   0
#define	LIMITER_MINMOD                                 1
#define LIMITER_VANLEER                                2
#define LIMITER_VANALBADA                              3
#define LIMITER_BARTH_JESPERSEN                        4
#define LIMITER_VENKATAKRISHNAN                        5
#define LIMITER_VENKATAKRISHNAN_CORRECTED              6

/**********************************************************************
 * CFD -- Flux Function Types.                                        *
 **********************************************************************/
#define	FLUX_FUNCTION_ROE                                   1
#define	FLUX_FUNCTION_HLLE                                  2
#define	FLUX_FUNCTION_HLLC                                  3

/**********************************************************************
 * CFD -- Local time-stepping/preconditioning.                        *
 **********************************************************************/

#define GLOBAL_TIME_STEPPING                                0
#define SCALAR_LOCAL_TIME_STEPPING                          1

/********************************************************
 * CFD -- Solver Types                                  *
 ********************************************************/

#define EXPLICIT                        0
#define IMPLICIT                        1

/*************************************************************************
 *                                                                       *
 * CFD -- Maximum length of codes and values of input parameter strings. *
 * Ideally this will be used for all input parameter classes.            *
 *                                                                       *
 *************************************************************************/

#define	INPUT_PARAMETER_MAX_LENGTH    128

/**********************************************************************
 * CFD -- Inline functions.                                           *
 **********************************************************************/

// minmod(x,y)
inline float  minmod(float x, float y)
   { return sgn(x)*max(float(ZERO),min(float(fabs(x)),sgn(x)*y)); }
inline double minmod(const double &x, const double &y)
   { return sgn(x)*max(ZERO,min(fabs(x),sgn(x)*y)); }

// minmod(x,y,z)
inline float  minmod(float x, float y, float z)
   { return sgn(x)*max(float(ZERO),min(float(fabs(x)),min(sgn(x)*y,sgn(x)*z))); }
inline double minmod(const double &x, const double &y, const double &z)
   { return sgn(x)*max(ZERO,min(fabs(x),min(sgn(x)*y,sgn(x)*z))); }

// minmod(x,y,z,w)
inline float  minmod(float x, float y, float z, float w)
   { return sgn(x)*max(float(ZERO),min(float(fabs(x)),min(sgn(x)*y,min(sgn(x)*z,sgn(x)*w)))); }
inline double minmod(const double &x, const double &y,
		     const double &z, const double &w)
   { return sgn(x)*max(ZERO,min(fabs(x),min(sgn(x)*y,min(sgn(x)*z,sgn(x)*w)))); }

// superbee(x,y)
inline float  superbee(float x, float y)
   { return sgn(x)*max(float(ZERO),max(min(float(TWO*fabs(x)),sgn(x)*y),
			           min(float(fabs(x)),float(TWO)*sgn(x)*y))); }
inline double superbee(const double &x, const double &y)
   { return sgn(x)*max(ZERO,max(min(TWO*fabs(x),sgn(x)*y),
			        min(fabs(x),TWO*sgn(x)*y))); }

// philimiter(x,y)
inline float  philimiter(float x, float y, float phi)
   { return sgn(x)*max(fabs(minmod(phi*x, y)), fabs(minmod(x, phi*y))); }
inline double philimiter(const double &x, const double &y,
			 const double &phi)
   { return sgn(x)*max(fabs(minmod(phi*x, y)), fabs(minmod(x, phi*y))); }

// vanleer(x,y)
inline float  vanleer(float x, float y)
   { return (fabs(x*y)+x*y)/(x+y+sgn(x+y)*sqr(TOLER)); }
inline double vanleer(const double &x, const double &y)
   { return (fabs(x*y)+x*y)/(x+y+sgn(x+y)*sqr(TOLER)); }

// vanalbada(x,y)
inline float  vanalbada(float x, float y, float epsi)
   { return (x*y+sqr(epsi))*(x+y)/(sqr(x)+sqr(y)+TWO*sqr(epsi)); }
inline double vanalbada(const double &x, const double &y,
		        const double &epsi)
   { return (x*y+sqr(epsi))*(x+y)/(sqr(x)+sqr(y)+TWO*sqr(epsi)); }

// Inline functions for AUSMplusUP flux calculation.
// M+1
inline double Mplus_1(double M)
  { return 0.5*(M + fabs(M)); }

// M-1
inline double Mminus_1(double M)
  { return 0.5*(M - fabs(M)); }

// M+2
inline double Mplus_2(double M)
  { return 0.25*sqr(M + 1.0); }

// M-2
inline double Mminus_2(double M)
  { return -0.25*sqr(M - 1.0); }

/********************************************************
 * CFD -- Define CFD structures and classes.            *
 ********************************************************/

/********************************************************
 * Class: CPUTime (CPU time for a process)              *
 *                                                      *
 * Member functions                                     *
 *      cpu_time0 -- Return CPU time used before last   *
 *                   update as a clock_t.               *
 *      cpu_time1 -- Return CPU time used at last       *
 *                   update as a clock_t.               *
 *      cput      -- Return total cput used by process. *
 *      reset     -- Resets CPU time counters.          *
 *      zero      -- Sets the CPU time to zero.         *
 *      update    -- Updates total CPU time for the     *
 *                   process.                           *
 *      sec       -- Returns total CPU time in seconds. *
 *      min       -- Returns total CPU time in minutes. *
 *      hrs       -- Returns total CPU time in hours.   *
 *                                                      *
 * Member operators                                     *
 *      T -- cpu time                                   *
 *                                                      *
 * T = T;                                               *
 * T = T + T;                                           *
 * T = T - T;                                           *
 * T = +T;                                              *
 * T = -T;                                              *
 * T += T;                                              *
 * T -= T;                                              *
 * T == T;                                              *
 * T != T;                                              *
 * cout << T; (output function)                         *
 * cin  >> T; (input function)                          *
 *                                                      *
 ********************************************************/
class CPUTime{
  private:
  public:
    clock_t cpu_time0,  // Track CPU time used so far as a clock_t.
            cpu_time1;  // To get the number of seconds used, 
                        // divide by CLOCKS_PER_SEC.
    double  cput;       // CPU time used by current process.
                        // Made public so can access them.

    /* Creation, copy, and assignment constructors. */
    CPUTime(void) {
       cpu_time0 = clock(); cpu_time1 = cpu_time0; 
       cput = ZERO;
    }

    CPUTime(const CPUTime &T) {
       cpu_time0 = T.cpu_time0; cpu_time1 = T.cpu_time1; 
       cput = T.cput; 
    }

    CPUTime(const double &cpu,
	    const clock_t &cpu0,
	    const clock_t &cpu1) {
       cpu_time0 = cpu0; cpu_time1 = cpu1; cput = cpu; 
    }

    /* Destructor. */
    // ~CPUTime(void);
    // Use automatically generated destructor.

    /* Reset CPU time counters. */
    void reset(void);

    /* Zero CPU time. */
    void zero(void);

    /* Update CPU time. */
    void update(void);

    /* CPU time in seconds. */
    double sec(void);
    double sec(void) const;

    /* CPU time in minutes. */
    double min(void);
    double min(void) const;

    /* CPU time in hours. */
    double hrs(void);
    double hrs(void) const;

    /* Assignment operator. */
    // CPUTime operator = (const CPUTime &T);
    // Use automatically generated assignment operator.

    /* Binary arithmetic operators. */
    friend CPUTime operator +(const CPUTime &T1, const CPUTime &T2);
    friend CPUTime operator -(const CPUTime &T1, const CPUTime &T2);
    
    /* Unary arithmetic operators. */
    friend CPUTime operator +(const CPUTime &T);
    friend CPUTime operator -(const CPUTime &T);

    /* Shortcut arithmetic operators. */
    friend CPUTime &operator +=(CPUTime &T1, const CPUTime &T2);
    friend CPUTime &operator -=(CPUTime &T1, const CPUTime &T2);
    
    /* Relational operators. */
    friend int operator ==(const CPUTime &T1, const CPUTime &T2);
    friend int operator !=(const CPUTime &T1, const CPUTime &T2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const CPUTime &T);
    friend istream &operator >> (istream &in_file,  CPUTime &T);

};

/********************************************************
 * CPUTime::reset -- Reset CPU time counters.           *
 ********************************************************/
inline void CPUTime::reset(void) {
    cpu_time0 = clock(); 
    cpu_time1 = cpu_time0; 
}

/********************************************************
 * CPUTime::zero -- Set total CPU time to zero.         *
 ********************************************************/
inline void CPUTime::zero(void) {
    cpu_time0 = clock(); 
    cpu_time1 = cpu_time0; 
    cput = ZERO;
}

/********************************************************
 * CPUTime::update -- Update total CPU time.            *
 ********************************************************/
inline void CPUTime::update(void) {
    cpu_time0 = cpu_time1;
    cpu_time1 = clock();
    cput = cput + double(cpu_time1-cpu_time0)/double(CLOCKS_PER_SEC);
}

/********************************************************
 * CPUTime::sec -- CPU time in seconds.                 *
 ********************************************************/
inline double CPUTime::sec(void) {
    return (cput);
}

inline double CPUTime::sec(void) const {
    return (cput);
}

/********************************************************
 * CPUTime::min -- CPU time in minutes.                 *
 ********************************************************/
inline double CPUTime::min(void) {
    return (cput/SIXTY);
}

inline double CPUTime::min(void) const {
    return (cput/SIXTY);
}

/********************************************************
 * CPUTime::hrs -- CPU time in hours.                   *
 ********************************************************/
inline double CPUTime::hrs(void) {
    return (cput/sqr(SIXTY));
}

inline double CPUTime::hrs(void) const {
    return (cput/sqr(SIXTY));
}

/********************************************************
 * CPUTime -- Binary arithmetic operators.              *
 ********************************************************/
inline CPUTime operator +(const CPUTime &T1, const CPUTime &T2) {
  return (CPUTime(T1.cput+T2.cput,T1.cpu_time0,T1.cpu_time1));
}

inline CPUTime operator -(const CPUTime &T1, const CPUTime &T2) {
  return (CPUTime(T1.cput-T2.cput,T1.cpu_time0,T1.cpu_time1));
}

/********************************************************
 * CPUTime -- Unary arithmetic operators.               *
 ********************************************************/
inline CPUTime operator +(const CPUTime &T) {
  return (CPUTime(T.cput,T.cpu_time0,T.cpu_time1));
}

inline CPUTime operator -(const CPUTime &T) {
  return (CPUTime(-T.cput,T.cpu_time0,T.cpu_time1));
}

/********************************************************
 * CPUTime -- Shortcut arithmetic operators.            *
 ********************************************************/
inline CPUTime &operator +=(CPUTime &T1, const CPUTime &T2) {
  T1.cput += T2.cput;
  T1.cpu_time0 = T1.cpu_time0;
  T1.cpu_time1 = T1.cpu_time1;
  return (T1);
}

inline CPUTime &operator -=(CPUTime &T1, const CPUTime &T2) {
  T1.cput -= T2.cput;
  T1.cpu_time0 = T1.cpu_time0;
  T1.cpu_time1 = T1.cpu_time1;
  return (T1);
}

/********************************************************
 * CPUTime -- Relational operators.                     *
 ********************************************************/
inline int operator ==(const CPUTime &T1, const CPUTime &T2) {
  return (T1.cput == T2.cput);
}

inline int operator !=(const CPUTime &T1, const CPUTime &T2) {
  return (T1.cput != T2.cput);
}

/********************************************************
 * CPUTime -- Input-output operators.                   *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const CPUTime &T) {
  out_file.setf(ios::scientific);
  out_file << " " << T.cput;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, CPUTime &T) {
  in_file.setf(ios::skipws);
  in_file >> T.cput;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * CFD -- External subroutines.                         *
 ********************************************************/
extern void Output_Progress_L2norm(const int Number_of_Time_Steps,
				   const double &Time,
				   const CPUTime &CPU_Time,
				   const double &Residual_L2_Norm,
				   const int First_Step,
				   const int Frequency);

extern void Output_Progress_L2norm(const int Number_of_Time_Steps,
				   const double &Time,
				   const CPUTime &CPU_Time,
				   const double &Residual_L2_Norm,
				   const int First_Step,
				   const int Frequency,
				   const int progress_character);

extern void Output_Progress_L2norm(const int Number_of_Time_Steps,
				   const double &Time,
				   const CPUTime &CPU_Time,
				   const double &Residual_L2_Norm,
				   const double &Ratio_Residual_L2_Norm,
				   const int First_Step,
				   const int Frequency);

extern void Output_Progress_L2norm(const int Number_of_Time_Steps,
				   const double &Time,
				   const CPUTime &CPU_Time,
				   const double &Residual_L2_Norm,
				   const double &Ratio_Residual_L2_Norm,
				   const int First_Step,
				   const int Frequency,
				   const int progress_character);

extern double Limiter_BarthJespersen(double *uQuad,
                                     const double &u0,
                                     const double &u0Min,
	      	                     const double &u0Max,
			             const int nQuad);

extern double Limiter_Venkatakrishnan(double *uQuad,
                                      const double &u0,
                                      const double &u0Min,
	      	                      const double &u0Max,
			              const int nQuad);

extern double Limiter_Venkatakrishnan_Modified(double *uQuad,
					       const double &u0,
					       const double &u0Min,
					       const double &u0Max,
					       const int nQuad);

extern double Limiter_VanLeer(double *uQuad,
                              const double &u0,
                              const double &u0Min,
	      	              const double &u0Max,
			      const int nQuad);

extern double Limiter_Minmod(double *uQuad,
			     const double &u0,
			     const double &u0Min,
			     const double &u0Max,
			     const int nQuad);

extern double Limiter_VanAlbada(double *uQuad,
                                const double &u0,
                                const double &u0Min,
	      	                const double &u0Max,
			        const int nQuad);

#endif /* _CFD_INCLUDED  */
