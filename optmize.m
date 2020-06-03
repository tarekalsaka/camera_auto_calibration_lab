clc
clear all;
format long g;
load('data.mat');
%=================================================================
% initial geuss  value for the intrinsic parameter
x0 = [A(1,1) A(1,2) A(1,3) A(2,2) A(2,3)];
%=================================================================
%  optimiser options
%  'TolX' (the step tolerance). defult 1e-6
%  'TolFun' (FunctionTolerance) ,defult 1e-6
%==================================================================
%=======minimize the cost function Mendonca-Cipolla Method CostFunctionMC
options = optimoptions(@lsqnonlin,'Algorithm',...
    'levenberg-marquardt','TolX',1e-7 );
OprimK_MC = lsqnonlin(@(x0) CostFunctionMC(Fs,x0),x0,[],[],options); 
finalMC  = [OprimK_MC(1) OprimK_MC(2) OprimK_MC(3);
    0 OprimK_MC(4) OprimK_MC(5); 0 0 1]
%==================================================================
%minimize the costfunction 
%The classical Kruppa’s equations costFunctionclassicK
options = optimoptions(@lsqnonlin,'Algorithm',...
    'levenberg-marquardt','TolFun',1e-8,'Tolx',1e-10);
OprimK_CK = lsqnonlin(@(x0) costFunctionclassicK(Fs,x0),x0,[],[],options);
finalCK  = [OprimK_CK(1) OprimK_CK(2) OprimK_CK(3)
    ; 0 OprimK_CK(4) OprimK_CK(5); 0 0 1]
%=========================================================================
%======minimize the cost function
%The Simplified Kruppa’sEquations costSimpifiedK
options = optimoptions(@lsqnonlin,'Algorithm',...
    'levenberg-marquardt','TolFun',1e-8,'Tolx',1e-20);
OprimK_SK = lsqnonlin(@(x0) costSimpifiedK(Fs,x0),x0,[],[],options); 
finalSK  = [OprimK_SK(1) OprimK_SK(2) OprimK_SK(3);
    0 OprimK_SK(4) OprimK_SK(5); 0 0 1]
%=====================================================================
%=======Calibration using the absolute dual Quadric
%find the Image Dual Abloute Quadric (IDAQ)
Q = IDAQlinearsolve();
Q0 = [Q(1,1) Q(1,2) Q(1,3) Q(1,4) Q(2,2) Q(2,3)...
    Q(2,4) Q(3,3) Q(3,4) Q(4,4)];
%minimizing the cost function Iterative methods %costFunItrativeMethod()
options = optimoptions(@lsqnonlin,'Algorithm',...
    'levenberg-marquardt','TolFun',1e-1,'Tolx',1e-9);
OprimK_IT= lsqnonlin(@(x0) costFunItrativeMethod(),x0,[],[],options);
finalIT = [OprimK_IT(1) OprimK_IT(2) OprimK_IT(3);
    0 OprimK_IT(4) OprimK_IT(5); 0 0 1]