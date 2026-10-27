# Gnuplot script for Lab 6 results
set terminal pngcairo enhanced font 'Arial,11' size 900,600
set datafile separator comma

# -------------------------------------------------------------------------------------
# PLOT 1: Convergence of Monte Carlo methods for \int_0^1 exp(x) dx
# -------------------------------------------------------------------------------------
set output 'fig_convergence_exp.png'
set title "Lab 6: Monte Carlo Convergence for \\int_0^1 e^x dx\nError vs Number of Samples (Log-Log Scale)" font ',13'
set xlabel "Sample Size N" font ',11'
set ylabel "Absolute Error |I_{est} - I_{exact}|" font ',11'
set grid mxtics mytics xtics ytics lc rgb "#d0d0d0" lt 1, lc rgb "#f0f0f0" lt 1
set logscale xy
set format x "10^{%T}"
set format y "10^{%T}"
set key top right box opaque

# Reference slope O(N^-0.5)
f_ref(x) = 0.5 / sqrt(x)

plot \
    'table3_convergence.csv' using 1:2 with linespoints pt 7 ps 1.2 lw 2 lc rgb "#0055d4" title "Sample Mean (Standard MC)", \
    'table3_convergence.csv' using 1:3 with linespoints pt 9 ps 1.2 lw 2 lc rgb "#d45500" title "Antithetic Variates", \
    'table3_convergence.csv' using 1:4 with linespoints pt 5 ps 1.2 lw 2 lc rgb "#009944" title "Control Variate (Linear)", \
    'table3_convergence.csv' using 1:5 with linespoints pt 11 ps 1.2 lw 2 lc rgb "#8800cc" title "Stratified Sampling", \
    f_ref(x) with lines dt 2 lw 2.5 lc rgb "#555555" title "Theoretical slope O(N^{-0.5})"

# -------------------------------------------------------------------------------------
# PLOT 2: Comparison of Monte Carlo vs Lab 5 Quadratures for \int_0^{\sqrt{3}} \sqrt{3-x^2} dx
# -------------------------------------------------------------------------------------
set output 'fig_comparison_task2.png'
set title "Comparison: Monte Carlo vs Lab 5 Newton-Cotes Quadratures\n\\int_0^{\\sqrt{3}} \\sqrt{3 - x^2} dx  (Quarter Circle, Singularity at x = \\sqrt{3})" font ',13'
set xlabel "Function Evaluations N" font ',11'
set ylabel "Absolute Error |I - Exact|" font ',11'
set grid mxtics mytics xtics ytics lc rgb "#d0d0d0" lt 1, lc rgb "#f0f0f0" lt 1
set logscale xy
set format x "10^{%T}"
set format y "10^{%T}"
set key top right box opaque

plot \
    'table2_comparison_lab5.csv' using (stringcolumn(2) eq "Средние прямоугольники (N=100)" || stringcolumn(2) eq "Средние прямоугольники (N=10000)" ? $3 : 1/0):6 with linespoints pt 7 ps 1.3 lw 2 lc rgb "#0066cc" title "Midpoint Rectangles [O(h^{1.5})]", \
    'table2_comparison_lab5.csv' using (stringcolumn(2) eq "Формула трапеций (N=100)" || stringcolumn(2) eq "Формула трапеций (N=10000)" ? $3 : 1/0):6 with linespoints pt 9 ps 1.3 lw 2 lc rgb "#cc6600" title "Trapezoidal Rule [O(h^{1.5})]", \
    'table2_comparison_lab5.csv' using (stringcolumn(2) eq "Формула Симпсона (N=100)" || stringcolumn(2) eq "Формула Симпсона (N=10000)" ? $3 : 1/0):6 with linespoints pt 5 ps 1.3 lw 2 lc rgb "#cc0000" title "Simpson's Rule [Degraded O(h^{1.5})]", \
    'table2_comparison_lab5.csv' using (stringcolumn(2) eq "Выборочное среднее (N=10^6)" ? $3 : 1/0):6 with points pt 4 ps 1.8 lw 2.5 lc rgb "#228833" title "MC Sample Mean (N=10^6)", \
    'table2_comparison_lab5.csv' using (stringcolumn(2) eq "Антитетические переменные (N=10^7)" ? $3 : 1/0):6 with points pt 6 ps 1.8 lw 2.5 lc rgb "#aa3377" title "MC Antithetic (N=10^7)"
