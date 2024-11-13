# IFC-Derivative-Signal-Processing

This project aims to automatically process impedance flow cytometry (IFC) raw data using a derivative-based method for peak detection. The notch filtering method has also been provided as an alternative signal processing technique for comparison.

## Project Overview

Impedance flow cytometry (IFC) is a powerful tool for analyzing cells or particles based on impedance signals. This project uses a derivative-based signal processing method to detect positive and negative peaks in the IFC data, allowing for accurate characterization of cell or particle passage events. The method is designed to handle real-world data, which often includes noise from various sources. 

The provided MATLAB code processes data files, detects peaks, and calculates peak-to-peak times to evaluate cell or particle characteristics.

## Main Features

- **Derivative Method for Peak Detection**: A derivative-based method is used to identify positive and negative peaks in the raw data stream. This method is particularly effective in handling noisy signals and provides accurate timing information for each peak event.
- **Notch Filtering Method**: In addition to the derivative method, a notch filtering method is available, allowing comparison between different processing techniques.
- **Automated Processing of Multiple Files**: The code processes multiple `.mat` data files stored in a specified folder, extracts relevant data for each file, and saves the results to CSV files.
- **Data Visualization**: The code generates plots of both the original and processed signals with detected peaks, enabling visual verification of the results.

## Code Structure

### Main Script: `PeakDetectionDerivative.m`, `PeakDetectionNotch.m`

The main MATLAB script performs the following steps:

1. **Load Data**: Loads raw data from `.mat` files in a specified folder.
2. **Preprocessing**: Detrends and normalizes the data to improve signal quality (derivative method does not require this step).
3. **Threshold Setting**: Sets a noise threshold to filter out unwanted peaks.
4. **Peak Detection**: Applies the derivative-based method to detect positive and negative peaks.
5. **Time Conversion**: Converts detected peak indices to time in milliseconds.
6. **Noise Handling**: Ignores files with excessive noise and skips saving if the peak count exceeds a specified threshold.
7. **Result Saving**: Saves the peak-to-peak times and amplitudes to CSV files for each channel and file.
8. **Data Visualization**: Plots the original signal and processed data with detected peaks for visual inspection.
