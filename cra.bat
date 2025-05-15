@echo off
g++ -o analysis analysis.cpp -lgdi32 -lgdiplus -mwindows
echo Compilation successful. Plotting Results...
analysis.exe
pause