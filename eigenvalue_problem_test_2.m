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
function [] = eigenvalue_problem_test_2
%
clc; clear eigenvalue_problem_test_2
%
A = [ 2, -1, 1;
     -1, 2, 0;
      1, 0, 6]; % matrix
%
lambda = eig(A) % The characteristic polynomial: p(lambda) = det(A - lambda*I) => -lambda^3 + 10*lambda^2 - 26*lambda + 16 = 0
                % The eigenvalues of A are the roots of the polynomial p(lambda)   
%
lambda = [0.8972
          2.8536
          6.2491];


%%%
return
end
