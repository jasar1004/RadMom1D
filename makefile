RADMOM1D_DIR=.
OTHERDIR_GRID=./Grid
OTHERDIR_CFD=./CFD
OTHERDIR_MEDIUM1DSTATE=./Medium1DState

LFLAGS = -g -lstdc++ -Wall -Wextra -lstdc++

# Objects for Grid
OBJS_GRID = $(OTHERDIR_GRID)/Cell1D.o

# Objects for CFD
OBJS_CFD = $(OTHERDIR_CFD)/CFD.o

# Objects for Medium1DState
OBJS_MEDIUM1DSTATE = $(OTHERDIR_MEDIUM1DSTATE)/Medium1DState.o

# Objects for RadMom1D
OBJS_RADMOM1D = $(RADMOM1D_DIR)/RadMom1DSolvers.o $(RADMOM1D_DIR)/RadMom1D_Mesh.o $(RADMOM1D_DIR)/RadMom1DInput_First_Order.o $(RADMOM1D_DIR)/RadMom1DInput_Third_Order.o $(RADMOM1D_DIR)/RadMom1DState_First_Order.o $(RADMOM1D_DIR)/RadMom1DState_Third_Order.o
OBJS_RADMOM1D_EXE = $(RADMOM1D_DIR)/RadMom1D.o

# Create the executable for RadMom1D
radMom1D: $(OBJS_RADMOM1D_EXE) $(OBJS_CFD) $(OBJS_GRID) $(OBJS_MEDIUM1DSTATE) $(OBJS_RADMOM1D)
	  g++ $(LFLAGS) -o radMom1D $(OBJS_RADMOM1D_EXE) $(OBJS_RADMOM1D) $(OBJS_CFD) $(OBJS_GRID) $(OBJS_MEDIUM1DSTATE)

# Objects for OBJS_RADMOM1D
RadMom1D.o: $(RADMOM1D_DIR)/RadMom1D.cc
		g++ $(LFLAGS) -c -o $(RADMOM1D_DIR)/RadMom1D.o $(RADMOM1D_DIR)/RadMom1D.cc

RadMom1DSolvers.o: $(RADMOM1D_DIR)/RadMom1DInput_First_Order.cc
		g++ $(LFLAGS) -c -o $(RADMOM1D_DIR)/RadMom1DSolvers.o $(RADMOM1D_DIR)/RadMom1DSolvers.cc

RadMom1D_Mesh.o: $(RADMOM1D_DIR)/RadMom1DInput_First_Order.cc
		g++ $(LFLAGS) -c -o $(RADMOM1D_DIR)/RadMom1D_Mesh.o $(RADMOM1D_DIR)/RadMom1D_Mesh.cc

RadMom1DInput_First_Order.o: $(RADMOM1D_DIR)/RadMom1DInput_First_Order.cc
		g++ $(LFLAGS) -c -o $(RADMOM1D_DIR)/RadMom1DInput_First_Order.o $(RADMOM1D_DIR)/RadMom1DInput_First_Order.cc
		
RadMom1DInput_Third_Order.o: $(RADMOM1D_DIR)/RadMom1DInput_Third_Order.cc
		g++ $(LFLAGS) -c -o $(RADMOM1D_DIR)/RadMom1DInput_Third_Order.o $(RADMOM1D_DIR)/RadMom1DInput_Third_Order.cc
		
RadMom1DState_First_Order.o: $(RADMOM1D_DIR)/RadMom1DState_First_Order.cc
		g++ $(LFLAGS) -c -o $(RADMOM1D_DIR)/RadMom1DState_First_Order.o $(RADMOM1D_DIR)/RadMom1DState_First_Order.cc
		
RadMom1DState_Third_Order.o: $(RADMOM1D_DIR)/RadMom1DState_Third_Order.cc
		g++ $(LFLAGS) -c -o $(RADMOM1D_DIR)/RadMom1DState_Third_Order.o $(RADMOM1D_DIR)/RadMom1DState_Third_Order.cc


# Objects for OBJS_GRID
Cell1D.o: $(OTHERDIR_GRID)/RadMom1DInput_First_Order.cc
		g++ $(LFLAGS) -c -o $(OTHERDIR_GRID)/Cell1D.o $(OTHERDIR_GRID)/Cell1D.cc

# Objects for OBJS_GRID
CFD.o: $(OTHERDIR_GRID)/CFD.cc
		g++ $(LFLAGS) -c -o $(OTHERDIR_GRID)/CFD.o $(OTHERDIR_GRID)/CFD.cc

# Objects for OBJS_MEDIUM1DSTATE
Medium1DState.o: $(OTHERDIR_MEDIUM1DSTATE)/Medium1DState.cc
		g++ $(LFLAGS) -c -o $(OTHERDIR_MEDIUM1DSTATE)/Medium1DState.o $(OTHERDIR_MEDIUM1DSTATE)/Medium1DState.cc

# Clean all object files in current directory (RadMom1D) and all associated subdirectories
clean:
	rm -rf *.o
	rm -rf radMom1D
	rm -rf $(OTHERDIR_GRID)/*.o
	rm -rf $(OTHERDIR_CFD)/*.o
	rm -rf $(OTHERDIR_MEDIUM1DSTATE)/*.o
