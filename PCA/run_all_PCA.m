%this script calls run_PCA for all the necessary inputs

dbstop if error
normalize = 1; %Normalization and centering is performed by cell and channel

%the following comment out code is already completed

%Shortwave:
useSW = 1;                
useLW = 0;
lwHiRes = 0;
run_PCA(useSW, useLW, lwHiRes, normalize);

%Longwave, low-res:
useSW = 0;                
useLW = 1;
lwHiRes = 0;
run_PCA(useSW, useLW, lwHiRes, normalize);

%Longwave, high-res:
useSW = 0;                
useLW = 1;
lwHiRes = 1;
run_PCA(useSW, useLW, lwHiRes, normalize);

%Shortwave and longwave, low-resolution:
useSW = 1;                
useLW = 1;
lwHiRes = 0;
run_PCA(useSW, useLW, lwHiRes, normalize);

%Shortwave and longwave, high-resolution:
useSW = 1;                
useLW = 1;
lwHiRes = 1;
run_PCA(useSW, useLW, lwHiRes, normalize);