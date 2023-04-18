/*!\file RadMom1DSolvers.h
  \brief 1D RadMom Equation Solvers. */

#ifndef _RADMOM_1D_SOLVER_INCLUDED
#define _RADMOM_1D_SOLVER_INCLUDED

/* Include 1D RadMom solution header file. */

#ifndef _RADMOM_1D_SINGLE_BLOCK_INCLUDED
#include "RadMom1D_Single_Block.h"
#endif // _RADMOM_1D_SINGLE_BLOCK_INCLUDED

/****************************************************************//**
 * Set default values for the input solution parameters
 * and then read user specified input values from the
 * specified input parameter file.
 ********************************************************************/
template<class cState, class pState>
int ReadInputFile(const char *Input_File_Name_ptr,
                  RadMom1D_Input_Parameters<cState, pState> &IP,
                  int &command_flag)
{
   // declares
   int error_flag(0);

   // The primary MPI processor processes the input parameter file.
   if (!command_flag) {
      cout << "\n Reading RadMom1D input data file `"
	       << Input_File_Name_ptr << "'." << endl;
   }

   error_flag = IP.Process_Input_Control_Parameter_File(Input_File_Name_ptr,
                                                        command_flag);

   if (!error_flag) {
      cout << IP << "\n";
      cout.flush();
   }

  // return error flag
  return error_flag;

}

/****************************************************************//**
 * Determine the L1, L2, and max norms of the solution
 * residual and update the CPU time for all processesors
 ********************************************************************/
template<class cState, class pState>
void ComputeNorms(RadMom1D_UniformMesh<cState, pState> *Soln,
                  double &residual_l1_norm,
                  double &residual_l2_norm,
                  double &residual_max_norm) {

  // L1 norm for all processors.
  residual_l1_norm = L1_Norm_Residual(Soln);

  // L2 norm for all processors.
  residual_l2_norm = L2_Norm_Residual(Soln);

  // Max norm for all processors.
  residual_max_norm = Max_Norm_Residual(Soln);

}

/********************************************************
 * Routine: RadMom1DSolver                               *
 *                                                      *
 * Computes solutions to 1D RadMom equations.            *
 *                                                      *
 ********************************************************/

template<class cState, class pState>
int RadMom1DSolver(const char *Input_File_Name_ptr) {

  /********************************************************
   * Local variable declarations.                         *
   ********************************************************/

  // 1D CFD input variables and parameters:
  RadMom1D_Input_Parameters<cState, pState> Input_Parameters;

  // Output file name pointer:
  char *Output_File_Name_ptr;

  // Output file stream:
  ofstream Output_File;

  /* Solution variables. */

  RadMom1D_UniformMesh<cState, pState> *Soln_ptr = NULL;

 /* Other local solution variables. */

  int number_of_time_steps,
      command_flag, error_flag, line_number;

  double time, dtime;

  double residual_l1_norm, residual_l2_norm, residual_max_norm(0);

  //! cpu time variables.
  CPUTime total_cpu_time;
  total_cpu_time.zero();

  /********************************************************
   * Set default values for the input solution parameters *
   * and then read user specified input values from the   *
   * specified input file.                                *
   ********************************************************/
  error_flag = ReadInputFile<cState, pState>(Input_File_Name_ptr, Input_Parameters, command_flag);
  if (error_flag) {
    cout << "\n RadMom1D ERROR: Unable to set default inputs.\n"
	 << flush;
    return error_flag;
  } // endif

  /*********************************************************
   * Create mesh and allocate RadMom1D solution variables.  *
   *********************************************************/

  execute_new_calculation: ;

  /* Allocate memory for 1D RadMom equation solution on
     uniform mesh. */

  cout << "\n Creating memory for RadMom1D solution variables.";
  cout.flush();

  Soln_ptr=Allocate<cState, pState>(Soln_ptr, Input_Parameters);

  if (Soln_ptr == NULL){
    cout << "\n RadMom1DSolvers::Allocate() Error! Probably not enough memory!";
    cout.flush();
    exit(1);
  }

  /* Create uniform mesh. */

  cout << "\n Creating uniform mesh.";
  cout.flush();

  Grid<cState, pState>(Soln_ptr,
                       Input_Parameters.X_Min,
                       Input_Parameters.X_Max,
                       Input_Parameters.Number_of_Cells_Idir);

  /********************************************************
   * Initialize RadMom1D solution variables.               *
   ********************************************************/

  /* Set the initial time level. */

  time = ZERO;
  number_of_time_steps = 0;

  /* Initialize the conserved and primitive state
     solution variables. */
  cout << "\n Prescribing RadMom1D initial data.";
  cout.flush();

  ICs<cState, pState>(Soln_ptr, Input_Parameters);

  cout << "\n Prescribing RadMom1D initial boundary conditions data.";
  cout.flush();

  BCs<cState, pState>(Soln_ptr, Input_Parameters);

  /********************************************************
   * Solve conservation form of 1D RadMom equations for    *
   * specified IBVP or BVP on uniform mesh.               *
   ********************************************************/

  continue_existing_calculation: ;

  if ((!Input_Parameters.Time_Accurate &&
          Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
          number_of_time_steps <= Input_Parameters.Maximum_Number_of_Time_Steps) ||
          (Input_Parameters.Time_Accurate &&
          Input_Parameters.Time_Max >= time)) {

    cout << "\n\n Beginning RadMom1D computations.\n\n";
    cout.flush();

     while (1) {
         /* Determine local and global time steps. */
         dtime = CFL<cState, pState>(Soln_ptr);

         // NORMS: Determine the L1, L2, and max norms of the
          // solution residual.
          ComputeNorms<cState, pState>(Soln_ptr, residual_l1_norm, residual_l2_norm, residual_max_norm);

          //screen
          Output_Progress_L2norm(number_of_time_steps,
                                 time*THOUSAND,
                                 total_cpu_time,
                                 residual_l2_norm,
                                 false, //first_step
                                 50);

          // CALCULATION CHECK: Check to see if calculations are
          // complete and if so jump of out of this infinite loop.
          if (!Input_Parameters.Time_Accurate &&
              (number_of_time_steps >= Input_Parameters.Maximum_Number_of_Time_Steps ||
              (residual_max_norm<= Input_Parameters.Min_Residual_Level && number_of_time_steps>0)))
              break;
          if (Input_Parameters.Time_Accurate && time >= Input_Parameters.Time_Max)
              break;

          // // LIMITER FREEZE: Freeze limiters as necessary
          // if (!first_step && Input_Parameters.Freeze_Limiter &&
          //     residual_l2_norm <= Input_Parameters.Freeze_Limiter_Residual_Level) {
          //     Freeze_Limiters(Local_SolnBlk, List_of_Local_Solution_Blocks);
          // } // endif

          // Update CPU time used for the calculation so far.
          total_cpu_time.update();

          // CALCULATION CHECK: Check to see if calculations are
          // complete and if so jump of out of this infinite loop.
          if (!Input_Parameters.Time_Accurate &&
              (number_of_time_steps >= Input_Parameters.Maximum_Number_of_Time_Steps ||
              (residual_max_norm<= Input_Parameters.Min_Residual_Level && number_of_time_steps>0)))
              break;
          if (Input_Parameters.Time_Accurate && time >= Input_Parameters.Time_Max)
              break;

        // 2. Apply boundary conditions for stage.
        BCs<cState, pState>(Soln_ptr, Input_Parameters);

         /* Update solution for next time step. */
         error_flag = dUdt_explicitRadMom_upwind<cState, pState>(Soln_ptr,
                                                                 Input_Parameters,
                                                                 dtime);

         if (error_flag) {
           cout << "\n PDES++ ERROR: RadMom1D solution error.";
           cout << "\n\nPDES++: Execution terminated.\n";
           cout.flush();
           return (error_flag);
         } /* endif */

         /* Update time and time step counter. */
         number_of_time_steps = number_of_time_steps + 1;
         time = time + Input_Parameters.CFL_Number*dtime;

     } /* endwhile */


     cout << "\n RadMom1D: Radiation Solver Stats:";
     cout << "\n RadMom1D: Total CPU time       ====> " << total_cpu_time.min() << " min";
     // cout << "\n RadMom1D: Total clock time     ====> " << difftime(end_explicit,start_explicit)/60.0 << " min";
     cout << "\n RadMom1D: Number of Time Steps ====> " << number_of_time_steps;
     cout << "\n RadMom1D: L1 Norm Residual     ====> " << residual_l1_norm;
     cout << "\n RadMom1D: L2 Norm Residual     ====> " << residual_l2_norm;
     cout << "\n RadMom1D: Max Norm Residual     ====> " << residual_max_norm;

    if (residual_max_norm > Input_Parameters.Min_Residual_Level) {
      cout << "\n RadMom1D ERROR: Failed to converge, Max Norm = " << residual_max_norm << " > Tolerance = " << Input_Parameters.Min_Residual_Level << ".";
    }
  } /* endif */

  // Recompute BCs
  BCs<cState, pState>(Soln_ptr, Input_Parameters);

  /********************************************************
   * Write 1D RadMom solution to output file.              *
   ********************************************************/

  while (1) {
     Input_Parameters.Get_Next_Input_Control_Parameter();
     command_flag = Input_Parameters.Parse_Next_Input_Control_Parameter();
     line_number = Input_Parameters.Line_Number;
     if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for 1D RadMom equation solution.
         cout << "\n Deallocating RadMom1D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr);
         // Execute new calculation.
         cout << "\n\n Starting a new calculation.";
         cout << Input_Parameters << "\n";

         goto execute_new_calculation;

     } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for 1D RadMom equation solution.
         cout << "\n Deallocating RadMom1D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr);
         // Close input data file.
         cout << "\n\n Closing RADMOM1D input data file.";
         Input_Parameters.Close_Input_File();
         // Terminate calculation.
         return (0);

     } else if (command_flag == CONTINUE_CODE) {
         // Reset maximum time step counter.
         Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;
         // Continue existing calculation.
         cout << "\n\n Continuing existing calculation.";
         cout << Input_Parameters << "\n";

         goto continue_existing_calculation;

     } else if (command_flag == WRITE_OUTPUT_CODE) {
         // Output solution data.
         Output_File_Name_ptr = Input_Parameters.Output_File_Name;

         cout << "\n Writing RadMom1D solution to output data file `"
              << Output_File_Name_ptr << "'.";
         Output_File.open(Output_File_Name_ptr, ios::out);
         if (Output_File.fail()) {
            cout << "\n PDES++ ERROR: Unable to open RadMom1D output data file.";
            cout << "\n\nPDES++: Execution terminated.\n";

            return (1);
         } /* endif */

            Output_Tecplot<cState, pState>(Soln_ptr,
                                           number_of_time_steps,
                                           time,
                                           Output_File);

         Output_File.close();
     } else if (command_flag == INVALID_INPUT_CODE ||
                command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         cout << "\nRadMom1D ERROR: Error reading RadMom1D data at line #"
              << -line_number  << " of input data file.";
         cout << "\n\nPDES++: Execution terminated.\n";

         return (line_number);
     } /* endif */
  } /* endwhile */

  /********************************************************
   * End of all RadMom1DSolver computations and I/O.       *
   ********************************************************/

  return (0);

}

#endif
