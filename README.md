# Wavelet-Quantile-on-Quantile-Regression-WQQR-with-P-values
To run the Wavelet Quantile-on-Quantile Regression (WQQR) with P-values, first ensure all required R packages (quantreg, np, waveslim, WaveletComp, ggplot2, plotly, and writexl) are installed and loaded. 
The script begins by loading the dataset (DATA.xlsx) and defining dependent (CO2) and independent (ICT) variables, which are then decomposed into short-, medium-, and long-term components using wavelet multi-resolution analysis (MRA). 
The Quantile Regression (QR) is performed across various quantiles (0.05 to 0.95), extracting slope coefficients and p-values, which are saved in Excel files for reference. 
Next, the Quantile-on-Quantile Regression (QQR) function estimates the relationship across quantiles of both variables, capturing nonlinear and asymmetric dependencies. 
The results are visualized using heatmaps, 3D surface plots, and comparative line graphs to contrast traditional quantile regression with WQQR estimates.
