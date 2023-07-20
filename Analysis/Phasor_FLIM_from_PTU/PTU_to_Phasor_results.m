
%  This code converts PTU files from Leica SP8 Falcon to Phasor to 
%  analyze fluorescence lifetime data. It outputs lifetime images with 
%  their respective  Phasor and an Excel File with the lifetime data.
%
%  This code performs batch analysis of PTU located in a folder.
%
%  This code calls the following functions:
%  - phasor_code.m
%  - PTU_LineScanRead_Phasor.m
%  - PTU_Read.m
%  - PTU_Read_Head.m
%  - PTU_to_Prehist.m
%  
%  The code runs on an i7-8700k with 16gb of RAM PC with the following
%  licenses:
%  -----------------------------------------------------------------------------------------------------
%  MATLAB Version: 9.14.0.2254940 (R2023a) Update 2
%  Operating System: Microsoft Windows 10 Enterprise Version 10.0 (Build 19045)
%  Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
%  -----------------------------------------------------------------------------------------------------
%  MATLAB                                                Version 9.14        (R2023a)
%  Curve Fitting Toolbox                                 Version 3.9         (R2023a)
%  Image Processing Toolbox                              Version 11.7        (R2023a)
%  Parallel Computing Toolbox                            Version 7.8         (R2023a)
%  Signal Processing Toolbox                             Version 9.2         (R2023a)
%  Statistics and Machine Learning Toolbox               Version 12.5        (R2023a)
%  Symbolic Math Toolbox                                 Version 9.3         (R2023a)
%
%  Part of this code is modified from code kindly provided by Enderlain
%  Laboratory in Goettingen.
%  Phasor code is modified from code by Pirmin Hubert Lakner
%  Lakner PH et al, Sci. Rep. 2017
%  (https://github.com/PirminLakner/Phasor_FLIM)
%  
%  (c) Mariano Gonzalez Pisfil, Updated July 2023 (Comments)

% Initialize program
clc
clear
close all

% Choose the folder where your PTU files. Only PTU files will be selcted
folder_name = uigetdir('\\nas.ads.mwn.de\ra56caz\MWN-PC\Dokumente\');

folder_FLIM = [folder_name '\*.ptu'];
E = dir(folder_FLIM);

tic % Start timer

% Start of Batch analysis, PTU converted to decay and image data. Once the
% processing is done, date is saved as .mat files so processing does not
% need to be redone from PTU files.

for i = 1:size(E,1)
    PTU_to_Prehist([E(i).folder '\' E(i).name], 25e6, 6);

  % PTU_to_Prehist(name, P, C)
  % name: name of the files processed
  % P: Photons batch size for processing. If your PC is slow of with small
  % RAM, choose a lower number 25 millions (25e6) works fine in the system
  % mentioned above
  % C: Cores in your PC procesor (Physical cores, not Logical Threads)

end

% ----------------------------------------------------------------------

% Now from the .mat files, the Phasor for each image in calculated. A
% threshold was sometimes necessary for our analysis. It can be chosen
% during the analysis. If different threshold are not wanted, it can be
% written on the software the threshold is not.

% Rescans the folder to list the .mat files
folder_FLIM = [folder_name '\*.mat'];
F = dir(folder_FLIM);

% Phasor code
phasor_code([], F(1).folder)

% Prompt to annouce the end of the analysis
fprintf('End of analysis, number of files analyzed: %f \n',i)

toc % Stop timer

close all