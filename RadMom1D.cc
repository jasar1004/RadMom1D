/*******************************************************************
 *******************************************************************
 **********         RadMom1D - Built standalone            *********
 ******                                                       ******   
 ****                                                           ****
 ******             Joachim Andre Raymond Sarr                ******
 **********                                                *********
 *******************************************************************
 ******************************************************************* 

 This program is a standalone version of the RadMom1D Solver
 that allows to solve the radiative transfer equation using
 moment closure techniques.

 The program is built standalone.

 *******************************************************************/

/* Include the header files defining various variable types and data
   structures, classes, functions, operators, global variables,
   etc.. */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

// Include cfrte1D header files.
#ifndef _RADMOM_1D_SOLVER_INCLUDED
#include "RadMom1DSolvers.h"
#endif // _RADMOM_1D_SOLVER_INCLUDED

/* Begin RadMom1D program. */

int main(int num_arg, const char *arg_ptr[]) {

  /********************************************************  
   * VARIABLE DECLARATIONS                                *
   ********************************************************/

  // Command line flags:  
  int error_flag;

  // Other local integer variables:
  int i;

  // Title of code:
  const char *program_title_ptr =
     "RadMom1D: Stand alone Solver for the radiative transfer equations in 1 space dimension.";
  
  // Input file name:
  const char *Input_File_Name_ptr = "radmom1D.in";

  // Equation type indicator:
  const char *Closure_Type_ptr = "M1";
  char Closure_Type[256];

  // Input file stream:
  ifstream Input_File;

  /********************************************************  
   * LOCAL VARIABLE INITIALIZATIONS                       *
   ********************************************************/

  /* Set default equation type. */

  strcpy(Closure_Type, Closure_Type_ptr);

  /********************************************************
   * PARSE COMMAND LINE ARGUMENTS                         *
   ********************************************************/

  /* Initialize command line flags. */
  error_flag = 0;

  /* Parse and interpret command line arguments.  Note that there
   * are several different possible arguments which are:
     1) -f name  uses "name" as the input data file rather than
                 the standard input data file "radmom1D.in". */

  if (num_arg >= 2) {
     for (i = 1; i <= num_arg - 1; ++i) {
      if (strcmp(arg_ptr[i],"-f") == 0) {

      } else if (strcmp(arg_ptr[i-1], "-f") == 0) {
         Input_File_Name_ptr = arg_ptr[i];
      } else {
         error_flag = 1;
      } /* endif */
    } /* endfor */
  } /* endif */


  /* Display command line argument error message and
     terminate the program as required. */

  if (error_flag) {
    cout << "\n RadMom1D ERROR: Invalid command line argument.\n";
    return (error_flag);
  } /* endif */


  /******************************************************************
   * DISPLAY THE PROGRAM TITLE AND VERSION INFORMATION AS REGUIRED. *
   ******************************************************************/

  cout << "**********************************************************************************************************" << '\n';
  cout << '\n' << program_title_ptr << '\n';
  cout << "Built by Joachim A. R. Sarr " << "\n";
  cout << '\n' << "**********************************************************************************************************" << '\n';
  cout.flush();

  /***********************************************************  
   * PERFORM REQUIRED CALCULATIONS.                          *
   ***********************************************************/
  if (strcmp(Closure_Type, "M1") == 0 || strcmp(Closure_Type, "P1") == 0) {
     error_flag = RadMom1DSolver<RadMom1D_cState_First_Order,
                                       RadMom1D_pState_First_Order>(Input_File_Name_ptr);
 } else if (strcmp(Closure_Type, "P3") == 0 ) {
    error_flag = RadMom1DSolver<RadMom1D_cState_Third_Order,
                                    RadMom1D_pState_Third_Order>(Input_File_Name_ptr);
 } /* endif */
  
  if (error_flag) {
     return (error_flag);
  } /* endif */


  /********************************************************  
   * TERMINATE PROGRAM EXECUTION                          *
   ********************************************************/
  cout << "\n\nRadMom1D: Execution complete.\n";

  //Ending properly
  return (0);

/* End RadMom1D program. */

}
