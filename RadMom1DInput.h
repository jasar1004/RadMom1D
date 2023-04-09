/////////////////////////////////////////////////////////////////////
/// \file RadMom1DInput.h
///
/// \author J. A. R. Sarr
///
/// This class defines all the input parameters read from the 
/// standard input and the input parameter file.
/// 
/////////////////////////////////////////////////////////////////////

#ifndef _RADMOM1D_INPUT_INCLUDED
#define _RADMOM1D_INPUT_INCLUDED

/* Include cfrte1D header files */
#ifndef _RADMOM1D_STATE_INCLUDED
#include "RadMom1DState.h" // Include 1D RadMom solution state header file
#endif // _RADMOM1D_STATE_INCLUDED

/////////////////////////////////////////////////////////////////////
// DEFINES
/////////////////////////////////////////////////////////////////////

// Define the structures and classes.
#define	INPUT_PARAMETER_LENGTH_RADMOM1D    128

/*!
 * \class RadMom1D_Input_Parameters
 * 
 *
 * @brief Definition and manipulation of 1D RadMom input variables.
 *
 */

template<class cState, class pState>
class RadMom1D_Input_Parameters;

template<class cState, class pState>
ostream &operator << (ostream &out_file, const RadMom1D_Input_Parameters<cState, pState> &IP);

template<class cState, class pState>  
istream &operator >> (istream &in_file, RadMom1D_Input_Parameters<cState, pState> &IP);

template<class cState, class pState>
class RadMom1D_Input_Parameters{
public:
    //! @name Input file parameters.
    //@{
    char     Input_File_Name[INPUT_PARAMETER_LENGTH_RADMOM1D]; //!< Input file name
    ifstream Input_File;                   //<! Input file stream:
    int      Line_Number;                  //<! Input file line number:
    //@}
    
    //@{ @name Time integration type indicator and related input parameters:
    char Time_Integration_Type[INPUT_PARAMETER_LENGTH_RADMOM1D]; //!< Time integrator string
    int i_Time_Integration;                                      //!< Time integrator flag
    int Time_Accurate;                                           //!< 0 false, 1 true
    int Local_Time_Stepping;                                     //!< time stepping flag
    int Maximum_Number_of_Time_Steps;                            //!< Max num explicit time steps
    double CFL_Number;                                           //!< CFL criteria
    double Time_Max;                                             //!< Max physical time
    
    //@{ @name Reconstruction parameters:
    char Reconstruction_Type[INPUT_PARAMETER_LENGTH_RADMOM1D]; //!< Reconstruction type string
    int i_Reconstruction;		                               //!< Reconstruction type flag
    char Limiter_Type[INPUT_PARAMETER_LENGTH_RADMOM1D];        //!< Limiter type string
    int i_Limiter;                                             //!< Limiter type flag
    int Freeze_Limiter;                                        //!< 0 false, 1 true
    double Freeze_Limiter_Residual_Level;
    double Min_Residual_Level;                                 //!< Residual cutoff
    static const int i_Residual_Variable = 1;                  //!< Residual variable to monitor
    static const int Number_of_Residual_Norms = 1;             //!< Number of residual norms
    //@}
    
    //@{ @name Moment closure type and related input parameters:
    char Moment_Closure_Type[INPUT_PARAMETER_LENGTH_RADMOM1D];  //!< String for M1, M2, P1 or P3
    int i_Moment_Closure;                                       //!< Flag for M1, M2, P1 or P3
    bool normalize;                                             //!< normalize solution
    //@}
    
    //! @name Flux function parameters
    //@{
    char Flux_Function_Type[INPUT_PARAMETER_LENGTH_RADMOM1D]; //!< Flux function string
    int i_Flux_Function;                                      //!< Flux function flag
    //@}
    
    //! @name Initial condition type indicator and related input parameters:
    char ICs_Type[INPUT_PARAMETER_LENGTH_RADMOM1D];       //!< IC type string
    int i_ICs;                                            //!< IC type flag
    char ICs_Intensity_or_Temperature[INPUT_PARAMETER_LENGTH_RADMOM1D]; // specify whether ICs should be based 
    // on given medium temperature or given intensity (string-based)
    int i_ICs_Flags; // specify whether ICs should be based on given medium temperature or given intensity (flag-based)
    pState Wo;                          //!< Reference primitive state
    cState Uo;                          //!< Reference conserved state
    Medium1D_State Mo;                  //!< Reference medium state
    double Pressure;                    //!< Pressure kPa
    double Temperature;                 //!< Temperature [K]
    double Reference_Temp;              //!< Reference temperature for scaling
    double Intensity;                   //!< Intensity
    double xco, xco2, xh2o, xo2, fsoot; //!< radiating gas composition
    //@}
    
    //@{ @name the corresponding case for the parallel plates tests:
    int Case;
    //@}
    
    //@{ @name Gas Parameters:
    double AbsorptionCoef;                              //!< Absorbsion coefficient [m^-1]
    double ScatteringCoef;                              //!< Scattering coefficient [m^-1]
    int i_ScatteringFunc;                               //!< Scattering phase function flag
    char ScatteringFunc[INPUT_PARAMETER_LENGTH_RADMOM1D]; //!< Scattering phase function string
    int i_AbsorptionModel;                              //!< Absorbsion model flag
    char AbsorptionModel[INPUT_PARAMETER_LENGTH_RADMOM1D]; //!< Absorption model string
    double OptThickness;
    double Albedo;
    //@}
    
    //@{ @name Boundary Conditions:
    double EastWallTemp, WestWallTemp;     //!< Wall temperature
    double EastWallEmiss, WestWallEmiss; //!< Wall emissivity
    //@}
    
    //! @name Grid Parameters
    int Number_of_Cells_Idir;                        //!< Number of cells in I-dir
    int Number_of_Ghost_Cells;                       //!< Number of ghost cells
    int Number_of_Blocks_Idir;                       //!< Number of blocks in I-dir
    double Box_Width;                                //!< Width of box

    double X_Min, X_Max; //!< coordinates for geometry's endpoints
    
    //@{ @name Boundary conditions:
    char Boundary_Conditions_Specified[INPUT_PARAMETER_LENGTH_RADMOM1D];
    int BCs_Specified;               //!< Flag indicating if boundary conditions are specified
    char Boundary_Conditions_Enforcement[INPUT_PARAMETER_LENGTH_RADMOM1D];
    int BCs_Enforcement;
    char BC_East_Type[INPUT_PARAMETER_LENGTH_RADMOM1D];   //!< Specified East BC string
    char BC_West_Type[INPUT_PARAMETER_LENGTH_RADMOM1D];   //!< Specified West BC string
    int BC_East;                                          //!< East BC flag
    int BC_West;                                          //!< West BC flag
    //@}
    
    //! @name Output Parameters
    //@{
    char Output_File_Name[INPUT_PARAMETER_LENGTH_RADMOM1D];          //!< Output file name
    char Grid_File_Name[INPUT_PARAMETER_LENGTH_RADMOM1D];            //!< Multi-block mesh definition input file names
    char Restart_File_Name[INPUT_PARAMETER_LENGTH_RADMOM1D];         //!< Restart file name
    char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_RADMOM1D];    //!< Next_Control_Parameter
    char Output_Format_Type[INPUT_PARAMETER_LENGTH_RADMOM1D];        //!< Output format type indicator
    int  i_Output_Format;                        //!< Flag indicating output format
    int  Restart_Solution_Save_Frequency;        //!< Restart solution frequency
    int  Output_Progress_Frequency;              //!< Output progress frequency
    //@}
    
    //! @name Parallel Domain Decomposition Parameters:
    //@{
    int Number_of_Processors;            //!< Total number of processors
    int Number_of_Blocks_Per_Processor;  //!< Number of blocks per processor
    int Number_of_Blocks;
    //@}
    
    //! Default constructor.
    RadMom1D_Input_Parameters(void);
    
    //! Default constructor.
    ~RadMom1D_Input_Parameters(void);

    void Setup_Radiative_Properties(void); 
    void SetupInputState(void);  //! Member function to setup conserved state Uo and medium state Mo
    void Open_Input_File(void);  //! open an input file
    void Close_Input_File(void); //! close an input file
    void Set_Default_Input_Parameters(void);
    // void Broadcast_Input_Parameters(void);        //! broadcast the input parameters
    // void Broadcast_Input_Parameters(MPI::Intracomm &Communicator, const int Source_CPU);
    void Get_Next_Input_Control_Parameter(void);  //! get the next input parameter
    int Parse_Next_Input_Control_Parameter(void); //! parse the next input parameter
    
    //! Output the name of the solver which this input parameters belong to.
    std::string Solver_Name(void){
        return "RadMom1D";
    }
    
    //! Process an input file
    int Process_Input_Control_Parameter_File(const char *Input_File_Name_ptr,
                                             int &Command_Flag);
    
    //@{ @name Input-output operators:
    friend ostream &operator << <cState, pState> (ostream &out_file, const RadMom1D_Input_Parameters<cState, pState> &IP);
    friend istream &operator >> <cState, pState> (istream &in_file, RadMom1D_Input_Parameters<cState, pState> &IP);
    //@}
};

/********************************************************
 * RadMom1D_Input_Parameters::RadMom1D_Input_Parameters()*
 * -->  Default Constructor                             *
 *******************************************************/
template<class cState, class pState>
inline RadMom1D_Input_Parameters<cState, pState>::RadMom1D_Input_Parameters(void){
    
}

/*********************************************************
 * RadMom1D_Input_Parameters::~RadMom1D_Input_Parameters()*
 * -->  Default Destructor                               *
 ********************************************************/
template<class cState, class pState>
inline RadMom1D_Input_Parameters<cState, pState>::~RadMom1D_Input_Parameters(void){
    pState::DeallocateStatic();
    cState::DeallocateStatic();
    Medium1D_State::DeallocateStatic();
}

/****************************************************************//**
* Output operator
********************************************************************/
template<class cState, class pState>
inline ostream &operator << (ostream &out_file,
                             const RadMom1D_Input_Parameters<cState, pState> &IP) {
    
    out_file << setprecision(6);
    
    out_file << endl 
            << endl << string(75,'*') 
            << endl << string(8,'*') << string(23,' ') << "RADMOM 1D INPUTS" << string(23,' ') << string(8,'*')
            << endl << string(75,'*');
            
    out_file << "\n\n Solving 1D Rte/RadMom equations (IBVP/BVP) on uniform mesh.";
    
    out_file << "\n - Input File Name = " << IP.Input_File_Name;
    
    if (IP.Time_Accurate) {
        out_file << "\n - Solution Type = UNSTEADY";
    } else {
        out_file << "\n - Solution Type = STEADY";
    }
    
    out_file << "\n - Time Integration = " << IP.Time_Integration_Type;
    
    out_file << "\n - Reconstruction = " << IP.Reconstruction_Type;
    
    out_file << "\n - Limiter = " << IP.Limiter_Type;   
    if (IP.Limiter_Type != LIMITER_ZERO && IP.Freeze_Limiter) {
        out_file << "\n     Freeze Limiter,  L2-norm < " 
                 << IP.Freeze_Limiter_Residual_Level;
    } 
    
    out_file << "\n - RADMOM Solver = " << IP.Moment_Closure_Type;
    
    if (IP.normalize) {
        out_file << "\n     Normalization = ON";
    } else {
        out_file << "\n     Normalization = OFF";
    }
    
    out_file << "\n - Absorption Model = " << IP.AbsorptionModel;
    out_file << "\n - Scattering Function = " << IP.ScatteringFunc;
    
    out_file <<"\n - Reference Temperature = " << IP.Reference_Temp;
    
    out_file << "\n - Initial Conditions = " << IP.ICs_Type;
    switch(IP.i_ICs) {
        case IC_CONSTANT :
        case IC_UNIFORM :
            out_file << "\n     Intensity = "  << IP.Intensity;
            out_file << "\n     Temperature (K) = " << IP.Temperature;
            out_file << "\n     Pressure (Pa) = " << IP.Pressure;
            out_file << "\n     Mixture = " 
                     << "xco = " << IP.xco << ",  "
                     << "xh2o = " << IP.xh2o << ",  "
                     << "xco2 = " << IP.xco2 << ",  "
                     << "xo2 = " << IP.xo2 << ",  "
                     << "fsoot = " << IP.fsoot;
            out_file << "\n     Absorption Cefficient (m^-1) = " << IP.AbsorptionCoef;
            out_file << "\n     Scattering Cefficient (m^-1) = " << IP.ScatteringCoef;
            break;
        default:
            out_file << "\n     Width (m) = " << IP.Box_Width;
            break;
    } // endswitch

    out_file << "\n     Width (m) = " << IP.Box_Width;
    
    if (IP.BCs_Specified) {
        out_file << "\n - Boundary conditions specified as: "
                 << "\n     BC_East = " << IP.BC_East_Type
                 << "\n     BC_West = " << IP.BC_West_Type;
    }

    out_file << "\n - Mesh Parameters";
    out_file << "\n     Number of Blocks I = " << IP.Number_of_Blocks_Idir;
    out_file << "\n     Number of Cells I = " << IP.Number_of_Cells_Idir;
    out_file << "\n     Number of Ghost Cells =  " << IP.Number_of_Ghost_Cells;
    
    out_file << "\n - Time Stepping Parameters";
    out_file << "\n     CFL Number = " << IP.CFL_Number;
    out_file << "\n     Maximum Time (ms) = " << IP.Time_Max*THOUSAND;
    out_file << "\n     Explicit Steps = " << IP.Maximum_Number_of_Time_Steps;
    
    out_file << "\n - CPU Parameters";
    out_file << "\n     Processors = " << IP.Number_of_Processors;
    out_file << "\n     Blocks Per Processor = " << IP.Number_of_Blocks_Per_Processor;
    
    out_file << "\n - Output Parameters";
    out_file << "\n     Output File Name = " << IP.Output_File_Name;
    out_file << "\n     Output Format = " << IP.Output_Format_Type;
    out_file << "\n     Restart Frequency = " << IP.Restart_Solution_Save_Frequency; 
    out_file << "\n     Output Progress Frequency = " << IP.Output_Progress_Frequency;
    
    out_file << endl << string(75,'*');
    
    return (out_file);
}

template<class cState, class pState>
inline istream &operator >> (istream &in_file,
                             RadMom1D_Input_Parameters<cState, pState> &IP) {
    return (in_file);
}

template<class cState, class pState>
void RadMom1D_Input_Parameters<cState,pState>::Setup_Radiative_Properties(void) {
    ScatteringCoef = Albedo * OptThickness;
    AbsorptionCoef = OptThickness - ScatteringCoef;
}

/******************************************************//**
 * Opens the appropriate input data file.             
 ********************************************************/
 template<class cState, class pState>
void RadMom1D_Input_Parameters<cState,pState>::Open_Input_File(void) {
    Input_File.open(Input_File_Name, ios::in);
    if (!Input_File.fail()) {
       Line_Number = 0;
       Input_File.setf(ios::skipws);
    } /* endif */

}

/******************************************************//**
 * Closes the appropriate input data file.      
 ********************************************************/
 template<class cState, class pState>
void RadMom1D_Input_Parameters<cState,pState>::Close_Input_File(void) {
    Input_File.unsetf(ios::skipws);
    Input_File.close();
}

/****************************************************************//**
 * Main function to setup the RadMom1D_State and Medium1D_State static
 * parameters.  This function is called everytime input parameters are
 * changed.
 ********************************************************************/
template<class cState, class pState>
void RadMom1D_Input_Parameters<cState,pState>::SetupInputState(void) 
{
  // deallocate
  Mo.Deallocate();
  Uo.Deallocate();
  Wo.Deallocate();
  
  // Setup static state variables
  Medium1D_State:: SetupStatic( i_AbsorptionModel, 
                                i_ScatteringFunc);
  
  pState::SetupStatic( i_ScatteringFunc, i_Moment_Closure, i_AbsorptionModel);
  
  cState::SetupStatic( i_ScatteringFunc, i_Moment_Closure, i_AbsorptionModel);

  // allocate
  Mo.Allocate();
  Uo.Allocate();
  Wo.Allocate();
  
  // initialize
  Mo.SetInitialValues( Temperature,
                       AbsorptionCoef,
                       ScatteringCoef);
  
  // Setup conserved and medium state
  if (i_ICs_Flags == RADMOM_ICS_INTENSITY) {
      Uo.Set_ICs_Intensity(Intensity);
  } else if (i_ICs_Flags == RADMOM_ICS_TEMPERATURE) {
      Uo.Set_ICs(Temperature);
  }
  
  Uo.closure_type = i_Moment_Closure;
  Wo = W(Uo);
  Wo.closure_type = i_Moment_Closure;
}

/******************************************************//**
 * Routine: Get_Next_Input_Control_Parameter            
 *                                                      
 * Get the next input control parameter from the input  
 * file.                                                
 *                                                      
 ********************************************************/
 template<class cState, class pState>
void RadMom1D_Input_Parameters<cState,pState>::Get_Next_Input_Control_Parameter(void) {
    int i;
    char buffer[256];
    
    Line_Number = Line_Number + 1;
    Input_File.getline(buffer, sizeof(buffer));
    
    if (Input_File.gcount() == sizeof(buffer)-1) {
        // if getline does not find a delimiter before size-1
        // characters then it sets the ifstream state to not
        // good.
        Input_File.clear(); 
        
        Input_File.ignore(10000, '\n');
        if (buffer[0] != '#') {
            cout << "\n***\n\nWarning: input file line " << Line_Number;
            cout << ": Line is more than " << sizeof(buffer) << " characters long. ";
            cout << "Ignoring rest of line.";
            cout << "\n\n***\n";
        }
    }
    
    i = 0;
    if (buffer[0] != '#') {
        while (1) {
            if (buffer[i] == ' ' || buffer[i] == '=' ) break;
            i = i + 1;
            if (i > strlen(buffer) ) break;
        } /* endwhile */
        buffer[i] = '\0';
    } /* endif */
    strcpy(Next_Control_Parameter, buffer);
}

/******************************************************//**
* Routine: Parse_Next_Input_Control_Parameter          
*                                                      
* Parses and executes the next input control parameter 
* from the input file.                                 
*                                                      
********************************************************/
template<class cState, class pState>
int RadMom1D_Input_Parameters<cState,pState>::Parse_Next_Input_Control_Parameter(void) {
    int i_command;
    char buffer[256];
    
    i_command = 0;
    
    if (strcmp(Next_Control_Parameter, "Time_Integration_Type") == 0) {
        i_command = 1;
        Get_Next_Input_Control_Parameter();
        strcpy(Time_Integration_Type, Next_Control_Parameter);
        if (strcmp(Time_Integration_Type, "Explicit_Euler") == 0) {
            i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
        } else {
            i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
        } /* endif */
    } else if (strcmp(Next_Control_Parameter, "Reconstruction_Type") == 0) {
        i_command = 2;
        Get_Next_Input_Control_Parameter();
        strcpy(Reconstruction_Type, Next_Control_Parameter);
        if (strcmp(Reconstruction_Type, "Green_Gauss") == 0) {
            i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
        } else if (strcmp(Reconstruction_Type, "Least_Squares") == 0) {
            i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
        } else {
            i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
        } /* endif */
    } else if (strcmp(Next_Control_Parameter, "Limiter_Type") == 0) {
        i_command = 3;
        Get_Next_Input_Control_Parameter();
        strcpy(Limiter_Type, Next_Control_Parameter);
        if (strcmp(Limiter_Type, "One") == 0) {
            i_Limiter = LIMITER_ONE;
        } else if (strcmp(Limiter_Type, "Zero") == 0) {
            i_Limiter = LIMITER_ZERO;
        } else if (strcmp(Limiter_Type, "VanLeer") == 0) {
            i_Limiter = LIMITER_VANLEER;
        } else if (strcmp(Limiter_Type, "VanAlbada") == 0) {
            i_Limiter = LIMITER_VANALBADA;
        } else if (strcmp(Limiter_Type, "Barth_Jespersen") == 0) {
            i_Limiter = LIMITER_BARTH_JESPERSEN;
        } else if (strcmp(Limiter_Type, "Venkatakrishnan") == 0) {
            i_Limiter = LIMITER_VENKATAKRISHNAN;
        } else {
            i_Limiter = LIMITER_VANLEER ;
        } /* endif */
    } else if (strcmp(Next_Control_Parameter, "Flux_Function_Type") == 0) {
        i_command = 4;
        Get_Next_Input_Control_Parameter();
        strcpy(Flux_Function_Type, Next_Control_Parameter);
        if (strcmp(Flux_Function_Type, "Roe") == 0) {
            i_Flux_Function = FLUX_FUNCTION_ROE;
        } else if (strcmp(Flux_Function_Type, "HLLE") == 0) {
            i_Flux_Function = FLUX_FUNCTION_HLLE;
        } else if (strcmp(Flux_Function_Type, "HLLC") == 0) {
            i_Flux_Function = FLUX_FUNCTION_HLLC;
        } else {
            i_Flux_Function = FLUX_FUNCTION_HLLE;
        } /* endif */
    } else if (strcmp(Next_Control_Parameter, "Moment_Closure_Type") == 0) {
        i_command = 5;
        Get_Next_Input_Control_Parameter();
        strcpy(Moment_Closure_Type, Next_Control_Parameter);
        if (strcmp(Moment_Closure_Type, "M1") == 0) {
            i_Moment_Closure  = MOMENT_CLOSURE_M1;
        } else if (strcmp(Moment_Closure_Type, "P1") == 0) {
            i_Moment_Closure  = MOMENT_CLOSURE_P1;
        } else if (strcmp(Moment_Closure_Type, "M2") == 0) {
            i_Moment_Closure  = MOMENT_CLOSURE_M2;
        } else if (strcmp(Moment_Closure_Type, "M2_Projection") == 0) {
            i_Moment_Closure  = MOMENT_CLOSURE_M2_PROJECTION;
        } else if (strcmp(Moment_Closure_Type, "P3") == 0) {
            i_Moment_Closure  = MOMENT_CLOSURE_P3;
        } else {
            i_Moment_Closure  = MOMENT_CLOSURE_M1;
        } /* endif */
    } else if (strcmp(Next_Control_Parameter, "ICs_Type") == 0) {
        i_command = 7;
        Get_Next_Input_Control_Parameter();
        strcpy(ICs_Type, Next_Control_Parameter);
        if (strcmp(ICs_Type, "None") == 0) {
            i_ICs = IC_NONE;
        } else if (strcmp(ICs_Type, "Uniform") == 0) {
            i_ICs = IC_UNIFORM;
        } else if (strcmp(ICs_Type, "Restart") == 0) {
            i_ICs = IC_RESTART;
        } else if (strcmp(ICs_Type, "Parallel_Plates") == 0) {
            i_ICs = IC_PARALLEL_PLATES;
        } else {
            std::cout << "\n ==> Unknown initial condition!";
            i_command = INVALID_INPUT_VALUE; exit(1);
        } /* endif */
    } else if (strcmp(Next_Control_Parameter, "ICs_Intensity_or_Temperature") == 0) {
        i_command = 8;
        Get_Next_Input_Control_Parameter();
        strcpy(ICs_Intensity_or_Temperature, Next_Control_Parameter);
        if (strcmp(ICs_Intensity_or_Temperature, "Intensity") == 0) {
            i_ICs_Flags = RADMOM_ICS_INTENSITY;
        } else if (strcmp(ICs_Intensity_or_Temperature, "Temperature") == 0) {
            i_ICs_Flags = RADMOM_ICS_TEMPERATURE;
        } else {
            i_ICs_Flags = RADMOM_ICS_TEMPERATURE;
        }
    } else if (strcmp(Next_Control_Parameter, "Output_File_Name") == 0) {
        i_command = 10;
        Get_Next_Input_Control_Parameter();
        strcpy(Output_File_Name, Next_Control_Parameter);
        strcat(Output_File_Name, ".dat");
        strcpy(Grid_File_Name, Next_Control_Parameter);
        strcat(Grid_File_Name, ".grid");
        strcpy(Restart_File_Name, Next_Control_Parameter);
        strcat(Restart_File_Name, ".soln");
    } else if (strcmp(Next_Control_Parameter, "Grid_File_Name") == 0) {
        i_command = 11;
        Get_Next_Input_Control_Parameter();
        strcpy(Grid_File_Name, Next_Control_Parameter);
        strcat(Grid_File_Name, ".grid");;
    } else if (strcmp(Next_Control_Parameter, "Restart_File_Name") == 0) {
        i_command = 12;
        Get_Next_Input_Control_Parameter();
        strcpy(Restart_File_Name, Next_Control_Parameter);
        strcat(Restart_File_Name, ".soln");
    } else if (strcmp(Next_Control_Parameter, "Number_of_Cells_Idir") == 0) {
        i_command = 13;
        Line_Number = Line_Number + 1;
        Input_File >> Number_of_Cells_Idir;
        Input_File.getline(buffer, sizeof(buffer));
        if (Number_of_Cells_Idir < 1) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "Number_of_Ghost_Cells") == 0) {
        i_command = 15;
        Line_Number = Line_Number + 1;
        Input_File >> Number_of_Ghost_Cells;
        Input_File.getline(buffer, sizeof(buffer));
        if (Number_of_Ghost_Cells < 1) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "Number_of_Blocks_Idir") == 0) {
        i_command = 16;
        Line_Number = Line_Number + 1;
        Input_File >> Number_of_Blocks_Idir;
        Input_File.getline(buffer, sizeof(buffer));
        if (Number_of_Blocks_Idir < 1) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "Number_of_Blocks") == 0) {
        i_command = 18;
        Line_Number = Line_Number + 1;
        Input_File >> Number_of_Blocks;
        Number_of_Blocks = Number_of_Blocks;
        Input_File.getline(buffer, sizeof(buffer));
        if (Number_of_Blocks < 1) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "Time_Accurate") == 0) {
        i_command = 19;
        Line_Number = Line_Number + 1;
        Input_File >> Time_Accurate;
        Input_File.getline(buffer, sizeof(buffer));
        if (Time_Accurate != 0 &&
            Time_Accurate != 1) {
            Time_Accurate = 0;
        }
        if (Time_Accurate) {
            Local_Time_Stepping = GLOBAL_TIME_STEPPING;
        } else {
            Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
        } /* endif */
    } else if (strcmp(Next_Control_Parameter, "Local_Time_Stepping") == 0) {
        i_command = 20;
        Line_Number = Line_Number + 1;
        Input_File >> Local_Time_Stepping;
        Input_File.getline(buffer, sizeof(buffer));
        if (Local_Time_Stepping != GLOBAL_TIME_STEPPING &&
            Local_Time_Stepping != SCALAR_LOCAL_TIME_STEPPING) {
            Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
        }
    } else if (strcmp(Next_Control_Parameter, "Maximum_Number_of_Time_Steps") == 0) {
        i_command = 21;
        Line_Number = Line_Number + 1;
        Input_File >> Maximum_Number_of_Time_Steps;
        Input_File.getline(buffer, sizeof(buffer));
        if (Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "CFL_Number") == 0) {
        i_command = 23;
        Line_Number = Line_Number + 1;
        Input_File >> CFL_Number;
        Input_File.getline(buffer, sizeof(buffer));
        if (CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "Box_Width") == 0) {
        i_command = 24;
        Line_Number = Line_Number + 1;
        Input_File >> Box_Width;
        Input_File.getline(buffer, sizeof(buffer));
        if (Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "OptThickness_Albedo") == 0) {
        i_command = 34;
        Line_Number = Line_Number + 1;
        Input_File >> OptThickness; if (OptThickness < ZERO) i_command = INVALID_INPUT_VALUE;
        Input_File >> Albedo; if (Albedo < ZERO) i_command = INVALID_INPUT_VALUE;
        Setup_Radiative_Properties();
        Input_File.getline(buffer, sizeof(buffer)); 
    } else if (strcmp(Next_Control_Parameter, "Absorption_Coefficient") == 0) {
        i_command = 35;
        Line_Number = Line_Number + 1;
        Input_File >> AbsorptionCoef;
        Input_File.getline(buffer, sizeof(buffer));
        if (AbsorptionCoef < ZERO) {
            i_command = INVALID_INPUT_VALUE;
        }  
    } else if (strcmp(Next_Control_Parameter, "Absorption_Model") == 0) {
        i_command = 36;
        Get_Next_Input_Control_Parameter();
        strcpy(AbsorptionModel, Next_Control_Parameter);
        if (strcmp(AbsorptionModel, "Gray") == 0) {
            i_AbsorptionModel = MEDIUM1D_ABSORB_GRAY;
        } else if (strcmp(AbsorptionModel, "SNBCK") == 0) {
            i_AbsorptionModel = MEDIUM1D_ABSORB_SNBCK;
        } else {
            i_command = INVALID_INPUT_VALUE;
        } /* endif */
    } else if (strcmp(Next_Control_Parameter, "Scattering_Coefficient") == 0) {
        i_command = 37;
        Line_Number = Line_Number + 1;
        Input_File >> ScatteringCoef;
        Input_File.getline(buffer, sizeof(buffer));
        if (ScatteringCoef < ZERO) {
            i_command = INVALID_INPUT_VALUE;
        }  
    } else if (strcmp(Next_Control_Parameter, "Scattering_Function") == 0) {
        i_command = 38;
        Get_Next_Input_Control_Parameter();
        strcpy(ScatteringFunc, Next_Control_Parameter);
        if (strcmp(ScatteringFunc, "Isotropic") == 0) {
            i_ScatteringFunc = RADIATION_SCATTER_ISO;
        } else if (strcmp(ScatteringFunc, "Linear") == 0) {
            i_ScatteringFunc = RADIATION_SCATTER_LINEAR;
        } else {
            i_command = INVALID_INPUT_VALUE;
        } /* endif */  
    } else if (strcmp(Next_Control_Parameter, "Gas_Temperature") == 0) {
        i_command = 39;
        Line_Number = Line_Number + 1;
        Input_File >> Temperature;
        Input_File.getline(buffer, sizeof(buffer));
        if (Temperature < ZERO) {
            i_command = INVALID_INPUT_VALUE;
      } 
    } else if (strcmp(Next_Control_Parameter, "Reference_Temperature") == 0) {
        i_command = 40;
        Line_Number = Line_Number + 1;
        Input_File >> Reference_Temp;
        Input_File.getline(buffer, sizeof(buffer));
        if (Reference_Temp <= ZERO) {
            i_command = INVALID_INPUT_VALUE;
        } 
    } else if (strcmp(Next_Control_Parameter, "Gas_Pressure") == 0) {
        i_command = 41;
        Line_Number = Line_Number + 1;
        Input_File >> Pressure;
        Input_File.getline(buffer, sizeof(buffer));
        if (Pressure < ZERO) {
            i_command = INVALID_INPUT_VALUE;
        }
    } else if (strcmp(Next_Control_Parameter, "Mixture") == 0) {
        i_command = 43;
        Line_Number = Line_Number + 1;
        Input_File >> xco;    if (xco   < ZERO) i_command = INVALID_INPUT_VALUE;
        Input_File >> xh2o;   if (xh2o  < ZERO) i_command = INVALID_INPUT_VALUE;
        Input_File >> xco2;   if (xco2  < ZERO) i_command = INVALID_INPUT_VALUE;
        Input_File >> xo2;    if (xo2   < ZERO) i_command = INVALID_INPUT_VALUE;
        Input_File >> fsoot;  if (fsoot < ZERO) i_command = INVALID_INPUT_VALUE;
        Input_File.getline(buffer, sizeof(buffer));
    } else if (strcmp(Next_Control_Parameter, "Intensity") == 0) {
        i_command = 44;
        Line_Number = Line_Number + 1;
        Input_File >> Intensity;
        Input_File.getline(buffer, sizeof(buffer));
        if (Intensity < ZERO) {
            i_command = INVALID_INPUT_VALUE;
        } 
    } else if (strcmp(Next_Control_Parameter, "Wall_Temperature") == 0) {
        i_command = 45;
        Line_Number = Line_Number + 1;
        Input_File >> EastWallTemp;  if (EastWallTemp < ZERO)  i_command = INVALID_INPUT_VALUE;
        Input_File >> WestWallTemp;  if (WestWallTemp < ZERO)  i_command = INVALID_INPUT_VALUE;
        Input_File.getline(buffer, sizeof(buffer));
    } else if (strcmp(Next_Control_Parameter, "Wall_Emissivity") == 0) {
        i_command = 46;
        Line_Number = Line_Number + 1;
        Input_File >> EastWallEmiss;  if (EastWallEmiss < ZERO)  i_command = INVALID_INPUT_VALUE;
        Input_File >> WestWallEmiss;  if (WestWallEmiss < ZERO)  i_command = INVALID_INPUT_VALUE;
        Input_File.getline(buffer, sizeof(buffer));
    } else if (strcmp(Next_Control_Parameter, "Time_Max") == 0) {
        i_command = 47;
        Line_Number = Line_Number + 1;
        Input_File >> Time_Max;
        Input_File.getline(buffer, sizeof(buffer));
        Time_Max = Time_Max/THOUSAND;
        if (Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "Number_of_Blocks_Per_Processor") == 0) {
        i_command = 48;
        Line_Number = Line_Number + 1;
        Input_File >> Number_of_Blocks_Per_Processor;
        Input_File.getline(buffer, sizeof(buffer));
        if (Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "Output_Format_Type") == 0) {
        i_command = 49;
        Get_Next_Input_Control_Parameter();
        strcpy(Output_Format_Type, Next_Control_Parameter);
        if (strcmp(Output_Format_Type, "Tecplot") == 0  ||
                   strcmp(Output_Format_Type, "tecplot") == 0  ||
                   strcmp(Output_Format_Type, "TECPLOT") == 0) {
            i_Output_Format = IO_TECPLOT;
        } else {
            i_Output_Format = IO_TECPLOT;
        } /* endif */
    } else if (strcmp(Next_Control_Parameter, "Case") == 0) {
        i_command = 54;
        Line_Number = Line_Number + 1;
        Input_File >> Case;
        Input_File.getline(buffer, sizeof(buffer));
        if (Case <= ZERO /*|| Case > 4*/) {
            i_command = INVALID_INPUT_VALUE;
        } /* endif */
    } else if (strcmp(Next_Control_Parameter, "Restart_Solution_Save_Frequency") == 0) {
        i_command = 57;
        Line_Number = Line_Number + 1;
        Input_File >> Restart_Solution_Save_Frequency;
        Input_File.getline(buffer, sizeof(buffer));
        if (Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter,"Output_Progress_Frequency") == 0) {
        i_command = 58;
        Line_Number = Line_Number + 1;
        Input_File >> Output_Progress_Frequency;
        Input_File.getline(buffer,sizeof(buffer));
        if (Output_Progress_Frequency < 1) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "Freeze_Limiter") == 0) {
        // Freeze_Limiter:
        i_command = 62;
        Line_Number = Line_Number + 1;
        Input_File >> Freeze_Limiter;
        Input_File.getline(buffer, sizeof(buffer));
        if (Freeze_Limiter < 0) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter, "Freeze_Limiter_Residual_Level") == 0) {
        // Freeze_Limiter_Residual_Level:
        i_command = 63;
        Line_Number = Line_Number + 1;
        Input_File >> Freeze_Limiter_Residual_Level;
        Input_File.getline(buffer, sizeof(buffer));
        if (Freeze_Limiter_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter,"Min_Residual_Level") == 0) {
        i_command = 82;
        Line_Number = Line_Number + 1;
        Input_File >> Min_Residual_Level;
        Input_File.getline(buffer,sizeof(buffer));
        if (Min_Residual_Level < 0) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(Next_Control_Parameter,"Boundary_Conditions_Specified") == 0) {
        i_command = 120;
        Get_Next_Input_Control_Parameter();
        strcpy(Boundary_Conditions_Specified,Next_Control_Parameter);
        if (strcmp(Boundary_Conditions_Specified,"ON") == 0) {
            BCs_Specified = ON;
        } else if (strcmp(Boundary_Conditions_Specified,"OFF") == 0) {
            BCs_Specified = OFF;
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    } else if (strcmp(Next_Control_Parameter,"Boundary_Conditions_Enforcement") == 0) {
        i_command = 121;
        Get_Next_Input_Control_Parameter();
        strcpy(Boundary_Conditions_Enforcement,Next_Control_Parameter);
        if (strcmp(Boundary_Conditions_Enforcement,"WEAK") == 0) {
            BCs_Enforcement = WEAK_BCS;
        } else if (strcmp(Boundary_Conditions_Enforcement,"STRONG") == 0) {
            BCs_Enforcement = STRONG_BCS;
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    } else if (strcmp(Next_Control_Parameter,"BC_East") == 0) {
        i_command = 126;
        Get_Next_Input_Control_Parameter();
        strcpy(BC_East_Type,Next_Control_Parameter);
        if (strcmp(BC_East_Type,"Reflection") == 0) {
            BC_East = BC_REFLECTION;
        } else if (strcmp(BC_East_Type,"Partial_Flux") == 0) {
            BC_East = BC_PARTIAL_FLUX;
        } else if (strcmp(BC_East_Type,"Partial_Moments") == 0) {
            BC_East = BC_PARTIAL_MOMENTS;
        } else if (strcmp(BC_East_Type,"Marshak") == 0) {
            BC_East = BC_MARSHAK;
        } else if (strcmp(BC_East_Type,"Fixed") == 0) {
            BC_East = BC_FIXED;
        } else if (strcmp(BC_East_Type,"Constant_Extrapolation") == 0) {
            BC_East = BC_CONSTANT_EXTRAPOLATION;
        } else if (strcmp(BC_East_Type,"Linear_Extrapolation") == 0) {
            BC_East = BC_LINEAR_EXTRAPOLATION;
        } else if (strcmp(BC_East_Type,"Characteristic") == 0) {
            BC_East = BC_CHARACTERISTIC;
        } else if (strcmp(BC_East_Type,"None") == 0) {
            BC_East = BC_NONE;
        } else if (strcmp(BC_East_Type,"Gray_Wall") == 0) {
            BC_East = BC_GRAY_WALL;
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    } else if (strcmp(Next_Control_Parameter,"BC_West") == 0) {
        i_command = 128;
        Get_Next_Input_Control_Parameter();
        strcpy(BC_West_Type,Next_Control_Parameter);
        if (strcmp(BC_West_Type,"Reflection") == 0) {
            BC_West = BC_REFLECTION;
        } else if (strcmp(BC_West_Type,"Partial_Flux") == 0) {
            BC_West = BC_PARTIAL_FLUX;
        } else if (strcmp(BC_West_Type,"Partial_Moments") == 0) {
            BC_West = BC_PARTIAL_MOMENTS;
        } else if (strcmp(BC_West_Type,"Marshak") == 0) {
            BC_West = BC_MARSHAK;
        } else if (strcmp(BC_West_Type,"Fixed") == 0) {
            BC_West = BC_FIXED;
        } else if (strcmp(BC_West_Type,"Constant_Extrapolation") == 0) {
            BC_West = BC_CONSTANT_EXTRAPOLATION;
        } else if (strcmp(BC_West_Type,"Linear_Extrapolation") == 0) {
            BC_West = BC_LINEAR_EXTRAPOLATION;
        } else if (strcmp(BC_West_Type,"Characteristic") == 0) {
            BC_West = BC_CHARACTERISTIC;
        } else if (strcmp(BC_West_Type,"None") == 0) {
            BC_West = BC_NONE;
        } else if (strcmp(BC_West_Type,"Gray_Wall") == 0) {
            BC_West = BC_GRAY_WALL;
        } else {
            i_command = INVALID_INPUT_VALUE;
        }
    } else if (strcmp(Next_Control_Parameter, "Execute") == 0) {
        i_command = EXECUTE_CODE;
    } else if (strcmp(Next_Control_Parameter, "Terminate") == 0) {
        i_command = TERMINATE_CODE;
    } else if (strcmp(Next_Control_Parameter, "Continue") == 0) {
        i_command = CONTINUE_CODE;
    } else if (strcmp(Next_Control_Parameter, "Write_Output") == 0) {
        i_command = WRITE_OUTPUT_CODE;
    } else if (strcmp(Next_Control_Parameter, "Write_Output_Cells") == 0) {
        i_command = WRITE_OUTPUT_CELLS_CODE;
    } else if (strcmp(Next_Control_Parameter,"Write_Output_Nodes") == 0) {
        i_command = WRITE_OUTPUT_NODES_CODE;
    } else if (strcmp(Next_Control_Parameter, "Write_Restart") == 0) {
        i_command = WRITE_RESTART_CODE;
    } else if (Next_Control_Parameter[0] == '#') {
        i_command = COMMENT_CODE;
    } else {
        i_command = INVALID_INPUT_CODE;
    } /* endif */
    
    // If it's still unknown then ignore it. 
    // This could be a bad idea if it was an unknown command 
    // as opposed to an unknown code.
    if (i_command == INVALID_INPUT_CODE) {
        cout << "\n***\n\nWarning: input file line " << Line_Number << ": ";
        cout << "ignoring unknown input code:\n";
        cout << "code: " << buffer;
        cout << "\nvalue: " << Next_Control_Parameter;
        cout << "\n\n***\n";
        i_command = COMMENT_CODE; // sure why not
    }
    
    if (!Input_File.good()) { i_command = INVALID_INPUT_VALUE; }
    
    /* Return the parser command type indicator. */
    return (i_command);
}

/******************************************************//**
 * Routine: Process_Input_Control_Parameter_File        
 *                                                      
 * Reads, parses, and executes the list of input        
 * control parameters from the standard input file.     
 *                                                      
 ********************************************************/
 template<class cState, class pState>
int RadMom1D_Input_Parameters<cState,pState>::Process_Input_Control_Parameter_File(const char *Input_File_Name_ptr,
                                                                                   int &Command_Flag) {
    
    int error_flag, line_number;
    
    /* Assign initial value for error indicator flag. */
    error_flag = 0;

    /* Assign default values to the input parameters. */
    Set_Default_Input_Parameters();
    
    /* Copy input file name (a string) to appropriate input parameter variable. */
    if (Input_File_Name_ptr != NULL) strcpy(Input_File_Name, Input_File_Name_ptr);
    
    /* Open the input file containing the input parameters. */
    Open_Input_File();
    error_flag = Input_File.fail();
    
    if (error_flag) {
        cout << "\n RadMom1D ERROR: Unable to open RadMom1D input input file.\n";
        return (error_flag);
    } /* endif */

    /* Read and parse control parameters contained in the input file. */
    while (1) {
        Get_Next_Input_Control_Parameter();
        Command_Flag = Parse_Next_Input_Control_Parameter();
        line_number = Line_Number;
        if (Command_Flag == EXECUTE_CODE) {
            // Setup conserved and medium state
            SetupInputState();
            break;
        } else if (Command_Flag == TERMINATE_CODE) {
            break;
        } else if (Command_Flag == INVALID_INPUT_CODE ||
                   Command_Flag == INVALID_INPUT_VALUE) {
            line_number = -line_number;
            cout << "\n RadMom1D ERROR: Error reading RadMom1D data at line #"
                 << -line_number  << " of input data file.\n";
            error_flag = line_number;
            break;
        } /* endif */
    } /* endwhile */
    
    // Recompute Number_of_Blocks_Per_Processor to ensure sufficient number of blocks are allocated for each processor
    Number_of_Blocks = Number_of_Blocks_Idir;
    Number_of_Blocks_Per_Processor = Number_of_Blocks/Number_of_Processors;
    

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */

    return (error_flag);

}

#endif /* _RADMOM1D_INPUT_INCLUDED  */
