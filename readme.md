# Image registration and data fusion 
This document aims to explain the general use of the program and its functions. 

## Goal of the program

## Use
1. Extract ROI for both LCM and HSI using getData. Need to be (roughly) same area. 
2. Register the two ROI's with each other using registration. 
3. Analyze the results with data_analysis



## Scripts
### Registration
Main script. Uses input of the extractRegion script to register the LCM and HSI datasets. 

### extractRegion
Extracts region of LCM and HSI data by user input 

### hsiContainer 
Class for convenicence to have all HSI data and functions in one place. 

### vk4data 
Class for convenicence to have all LCM data and functions in one place. 

### Theory_test
Compare results with *simple* theory. 


### data_analysis
Take registration outputs and analyze 4D dataset. 
Compare height with wavelengths. 