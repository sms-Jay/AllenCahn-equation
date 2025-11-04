@echo off
title Allen-Cahn Equation Compiler
color 0A

echo ========================================
echo    Allen-Cahn Equation Solver Compiler
echo ========================================
echo.

echo [1/2] Compiler
g++ -std=c++11 -O3 -fopenmp -I "include" "src/nontruncated_admm_parallel.cpp" -o "src/nontruncated_admm_parallel.exe"
if %errorlevel% equ 0 (
    echo [OK] ADMM compiled successsfully
) else (
    echo [ERROR] ADMM compiled NOT successsfully
    goto :error
)

echo.
echo [2/2] Compiler
g++ -std=c++11 -O3 -fopenmp -I "include" "src/newton_parallel.cpp" -o "src/newton_solver.exe"
if %errorlevel% equ 0 (
    echo [OK] Newton compiled successsfully
) else (
    echo [ERROR] Newton compiled NOT successsfully
    goto :error
)

echo.
echo ========================================
echo   All compiled!
echo ========================================
echo.
echo 
echo   src\nontruncated_admm_parallel.exe
echo   
echo   src\newton_solver.exe
echo.
pause
exit /b 0

:error
echo.
echo ========================================
echo   Something wrong!
echo ========================================
echo.
pause
exit /b 1