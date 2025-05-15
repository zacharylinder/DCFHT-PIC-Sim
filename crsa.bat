@echo off
g++ -funroll-all-loops -march=native -O3 -ffast-math -o simulation simulator.cpp
echo Compilation successful. Running simulation...
simulation.exe
echo.
echo Simulation successful. Compiling Analysis...
g++ -o analysis analysis.cpp -lgdi32 -lgdiplus -mwindows
echo Compilation successful. Plotting Results...
analysis.exe
pause
