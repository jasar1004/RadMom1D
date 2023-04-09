/*!\file RadMom1DInput.cc
  \brief Subroutines for 1D Radiative Moment Input Classes. */

/* Include RADMOM header files */ 
#include "RadMom1DInput.h" /* Include 1D RadMom input parameter header file. */

// static variables
template<>
const int RadMom1D_Input_Parameters<RadMom1D_cState_First_Order,
                                    RadMom1D_pState_First_Order>::i_Residual_Variable;

template<>
const int RadMom1D_Input_Parameters<RadMom1D_cState_First_Order,
                                    RadMom1D_pState_First_Order>::Number_of_Residual_Norms;
                                    
/******************************************************//**
 * Routine: Set_Default_Input_Parameters               
 *                                                     
 * Assigns default values to the input parameters.     
 *                                                     
 ********************************************************/
template<>
void RadMom1D_Input_Parameters<RadMom1D_cState_First_Order,RadMom1D_pState_First_Order>::Set_Default_Input_Parameters(void) {

    strcpy(Input_File_Name, "RadMom1D.in");

    strcpy(Time_Integration_Type, "Explicit_Euler");
    i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    Time_Accurate = 0;
    Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
    Maximum_Number_of_Time_Steps = 100;
    CFL_Number = 0.3;
    Time_Max = ZERO;

    // Reconstruction type:
    strcpy(Reconstruction_Type, "Least_Squares");
    i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
    
    // Limiter type:
    strcpy(Limiter_Type, "Barth_Jespersen");
    i_Limiter = LIMITER_BARTH_JESPERSEN;

    // Flux function:
    strcpy(Flux_Function_Type, "HLLE");
    i_Flux_Function = FLUX_FUNCTION_HLLE;

    // Moment function:
    strcpy(Moment_Closure_Type, "M1");
    i_Moment_Closure = MOMENT_CLOSURE_M1;
    normalize = false;

    // Initial conditions:
    strcpy(ICs_Type, "Uniform");
    strcpy(ICs_Intensity_or_Temperature, "Temperature");
    i_ICs = IC_UNIFORM;
    i_ICs_Flags = RADMOM_ICS_TEMPERATURE;
    Intensity      = ZERO;
    Temperature    = THOUSAND;       //[K]
    // Pressure       = PRESSURE_STDATM;//[Pa]
    xco            = 0.01;
    xh2o           = 0.2;
    xco2           = 0.1;
    xo2            = 0.0;
    fsoot          = 0.0;
    Reference_Temp = Temperature;
    Case = 0;
    
    // Gas parameters
    i_AbsorptionModel = MEDIUM1D_ABSORB_GRAY;
    strcpy(AbsorptionModel, "Gray");
    i_ScatteringFunc = RADIATION_SCATTER_ISO;
    strcpy(ScatteringFunc, "Isotropic");
    AbsorptionCoef = ONE;
    ScatteringCoef = ZERO;
    OptThickness = ONE;
    Albedo = ZERO;

    // boundary conditions
    EastWallTemp = ZERO;
    WestWallTemp = ZERO;
    EastWallEmiss = ZERO;      
    WestWallEmiss = ZERO;  

    // Grid parameters:
    Box_Width = ONE;
    X_Min = ZERO;
    X_Max = ONE;
    Number_of_Cells_Idir = 100;
    Number_of_Ghost_Cells = 2;
    Number_of_Blocks_Idir = 1;
    Number_of_Blocks = 1;

    // Boundary conditions:
    strcpy(Boundary_Conditions_Specified,"OFF");
    BCs_Specified = OFF;
    BC_East = BC_GRAY_WALL;
    BC_West = BC_GRAY_WALL;
    
    strcpy(Boundary_Conditions_Enforcement,"WEAK");
    BCs_Enforcement = WEAK_BCS;
    
    strcpy(BC_East_Type,"None");
    strcpy(BC_West_Type,"None");
    BC_East  = BC_NONE;
    BC_West  = BC_NONE; 

    // Default output file names and parameters:
    strcpy(Output_File_Name, "outputfile.dat");

    strcpy(Grid_File_Name, "gridfile.grid");
    
    strcpy(Restart_File_Name, "restartfile.soln");

    strcpy(Output_Format_Type, "Tecplot");
    i_Output_Format = IO_TECPLOT;
    Restart_Solution_Save_Frequency = 10000;
    Output_Progress_Frequency = 50; 

    // Input_file parameters:
    strcpy(Next_Control_Parameter, " ");

    Line_Number = 0;

    // Number_of_Blocks_Per_Processor = 10;

    // Freezing_Limiter
    Freeze_Limiter = 0;
    Freeze_Limiter_Residual_Level = 1e-4;
    Min_Residual_Level = TOLER;

    // Setup conserved and medium state
    SetupInputState();
}
