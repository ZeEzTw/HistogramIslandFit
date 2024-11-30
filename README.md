# HistogramIslandFit
Overview

This C++ code processes histograms stored in ROOT files and performs Gaussian fitting on regions of interest (ROIs) to find the centers of two distinct peaks. The program works with 32 detectors and processes data for 2 energy values. It uses ROOT for histogram manipulation and data visualization.
Key Functionality

    Region Search:
        For each pair of detectors, after applying a formula, it identifies two regions of interest (peaks).
        The algorithm searches within a defined area for the largest clusters of bins that correspond to the peaks.

    Gaussian Fitting:
        Applies a 2D Gaussian fit to the identified regions (peaks) to find their centers.
        Draws circles around the identified peak centers on the histograms.

    Line Fitting:
        Draws a line between the two identified peak centers and performs a linear fit to calculate the slope and intercept.

    Calculation:
        Calculates a delta value (delta_ij) based on the fit parameters and stores it in an output file.
        Also calculates the reverse delta (delta_ji) and saves the results.

    Output:
        Generates a ROOT file with histograms and fitted lines.
        Saves the fit parameters and calculated deltas in a text file.

Output Files

    output/results.txt: Contains the calculated delta values (delta_ij and delta_ji).
    output/histograms_with_lines_and_slopes.root: ROOT file containing the histograms with fitted Gaussian curves and lines.
    Individual ROOT files for each pair of detectors (histogram_with_line_i_j.root).

## How to Run It

1. Open ROOT with the command:
   
        root --web=off
 
2. Load and execute the code using:

       .L findBiggestBinContent.cpp
       
