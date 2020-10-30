# RIS-enabled-SISO-Localization

## Paper
This project is based on the paper "SISO RIS-Enabled Joint 3D Downlink Localizationand Synchronization". 

## Problem statement

We consider a wireless system with a single-antenna BS and UE, as well as an RIS.   The user receives the transmitted signal directly from  BS (LOS path) and also from  RIS (reflected path). We assume that the position of the BS and  RIS as well as the rotation matrix corresponding to the RIS coordinate system are known On the other hand, the UE position is unknown and should be estimated. In addition,  we consider an asynchronous scenario, where the user clock bias should also be estimated as a nuisance parameter.

## Parameters
The parameters are the same as those in the paper, except the number of iteration that is 100 (5000 in the paper). They can be adjusted inside the file GetConfig.m.

## Runtime 
Upon running the main.m file all the figures in the paper are generated in less than 2 hours . Remember that the number of iterations are much lower (100 vs 5000 in paper), which can make the result look different. You can adjust the number of iteration by changing "config.NoiseSampleNum".

## Matlab version
The code has been written in MATLAB 2018a
