%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes an eigenvalue of matrix A. The eigenvalues are
% obtained as roots of the characteristic polynomial: p(lambda) = det(A - lambda*I) = 0 [1].
%
% [1] Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 29, 2025
%
%%% 2025.12.29
function [] = eigenvalue_problem_test_1
%
clc; clear eigenvalue_problem_test_1
%
A = [18, 10;
     10, 13]; % matrix
%
lambda = eig(A) % The characteristic polynomial: p(lambda) = det(A - lambda*I) = lambda^2 - 31*lambda + 134 = 0 => lambda = (31+/- 5*sqrt(5))/2
                % The eigenvalues of A are the roots of the polynomial p(lambda)   
%
lambda = [5.1922
         25.8078];


%%%
return
end
