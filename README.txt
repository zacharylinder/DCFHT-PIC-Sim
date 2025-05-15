
COMPILE SCRIPTS

    simulator.cpp: 
        g++ -funroll-all-loops -march=native -O3 -ffast-math -o simulation simulator.cpp
    analysis.cpp: 
        g++ -o analysis analysis.cpp -lgdi32 -lgdiplus -mwindows

FILE STRUCTURE

    progressbar.h - Header file for handling the visual progress bar updating (doesn't need to be changed)
    simulator.cpp - Main C++ file for performing all calculations and data logging
    bField.fld - B Field data (from ANSYS sim)
    eField.fld - E Field data (from ANSYS sim)
    discrete.txt - The amount of particles each cell sees over the course of the simulation
        - Data to be logged for cell plot can be defined in settings.txt under "PLOT", which indicates 
          the number of equally spaced timestamps it logs
    log.txt - The positions of all particles seen over the course of the simulation
        - Data to be logged for particle plot can be defined in settings.txt under "POINTS", which indicates
          the number of equally spaced timestamps it logs.
    settings.txt - Settings to define all independent variables
    crs.bat - Compiles simulator.cpp, then runs simulation.exe, then runs analysis.exe
        - Will only work if you have the g++ compiler and the necessary C++ libraries as indicated in the 
          #include section of simulator.cpp
    cra.bat - Compiles analysis.cpp, then runs analysis.exe
        - Will only work if you have the g++ compiler and the necessary C++ libraries as indicated in the 
          #include section of simulator.cpp
    run.bat - Runs simulation.exe, then runs analysis.exe
    simulation.exe - A result of compiling simulator.cpp. Runs the simulation (computation & logging only)
    analysis.cpp - Displays all logged data in a positions/heatmap plot and a cell-based plot in W11 window
    icon.ico - Window icon at top-left of screen

SETTINGS

    To edit simulation settings, you can adjust the values inside settings.txt (do not change names)
        PARTICLES - Number of Particles to simulate
        DT - The time inbetween each timestep
        TIMESTEPS - Number of timesteps to simulate (simulation time would be DT * TIMESTEPS)
        POINTS - Number of positions per particle to record in particle plot
        QOM - Charge-To-Mass Ratio (this is a constant)
        RADIUS - Radius of the outlet
        INLET - Radius of the inlet 
        DIV - Divergence Angle
        PLOT - Number of timestamps to take for cell densities 

PSEDUOCODE

    0 - Import necessary data (field data, settings, etc.)
    1 - Initialize particles at top of outlet with initial position and velocities
    2 - Define geometrical parameters of thruster given initial conditions in settings.txt
    3 - Run for-loop to iterate particle parameters
        a - Determine EB field contribution by particles onto the grid's EB field
            i - Take a weighted vector by the 4 quadrants the particles makes with a cell
            ii - Multiply the particle's EB field by the appropriate weight and append it to the field
        b - Determine the interpolated EB field by the grid onto the particles' EB fields
            i - Take a weighted vector by the 4 quadrants the particles makes with a cell
            ii - Multiply the 4 grid EB fields by the appropriate weights and define it to the particle
            iii - Record cell density if necessary
        c - Perform Maxwells equations to achieve new position and velocity given the particles' EB fields
            i - Check for wall collisions, in which case the velocity vector is reflected across 
                the wall's normal (if so, reinitialize)
        d - Record particle positions if necessary, and updated progress bar
    4 - Post-processing of simulation results

