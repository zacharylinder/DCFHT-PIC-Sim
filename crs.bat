@echo off
g++ -funroll-all-loops -march=native -O3 -ffast-math -o simulation simulator.cpp
echo Compilation successful. Running simulation...
simulation.exe
echo.
echo Simulation successful. Plotting Results...
analysis.exe
