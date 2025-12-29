%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves the one-dimensional (dim) Schrodinger equation (TISE) for the
% 1D soft-Coulomb (V(x) = -1/sqrt(1+x^2)) [1], using finite difference method with Dirichlet boundary condition [2]. 
%
% The atomic units are used (\hbar = m = e = 1). 
% The (TISE): H(x)*psi(x) = En*psi(x) => (-0.5*d^2/dx^2 + V(x))*psi(x) = En*psi(x) 
% The Dirichlet boundary condition: psi(0) = psi(L) = 0.         
%
% [1] Q Su and JH Eberly. Model atom for multiphoton physics. Physical Review A, 44(9):5997, 1991.
% [2] Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 29, 2025
%
%%% 2025.12.29
function [] = one_dim_soft_Coulomb_time_indepedent_Schrodinger_eq
%
clc; clear one_dim_soft_Coulomb_time_indepedent_Schrodinger_eq
%
% grid and initial condition
a = -25.; % x(0)
b = 25.; % x(N+1)
N = 512;  % number of grid point of x axis
dx = (b-a)/N; % step size in x
% ---
x = zeros(N+1,1); % total number of points is N+1
for i = 1:N+1
    x(i) = a + (i-1)*dx;
end
% ---
lambda = 1./dx^2 ;
%  matrix equation is
u_mat = zeros(N+1,N+1);
for i = 2:N
    u_mat(i,i-1) = lambda;
    u_mat(i,i) = -2.*lambda;
    u_mat(i,i+1) = lambda;
end
u_mat(1,1) = -2.*lambda; u_mat(N+1,N+1) = -2.*lambda ; 
u_mat(1,2) = lambda; u_mat(N+1,N) = lambda;
%u_mat;
%%%
V_soft_Coulomb_pot = -1./sqrt(1.+x.^2);
%
H0_ham = -0.5.*u_mat + diag(V_soft_Coulomb_pot); % Hamiltonian (H0)
%
[Vec,En] = eig(H0_ham(2:N,2:N));                                     % Eigenvalue problem with the Dirichlet boundary condition
En = diag(En);
[foo, ij] = sort(En);
En = En(ij);
[En(1),En(2),En(3),En(4),En(5)]
%
ans = [-0.6699   -0.2750   -0.1515   -0.0927   -0.0632];
%   = [-0.6698   -0.2749   -0.1515   -0.0927   -0.0635 ]
%
   
%%%
%           
Vec = Vec(:,ij);                       % The unnormalized eigenfunctions
V1 = Vec(:,1);                         % The unnormalized eigenfunction for the ground state,
%V1 = [0.,;V1,;0.];
n_c = sum(V1.*V1.*dx);
V1 = 1./sqrt(n_c).*V1;
%
V2 = Vec(:,2);                         % The unnormalized eigenfunction for the 1st excited state,
%V1 = [0.,;V1,;0.];
n_c2 = sum(V2.*V2.*dx);
V2 = 1./sqrt(n_c2).*V2;
%
V3 = Vec(:,3);                         % The unnormalized eigenfunction for the 2nd excited state,
%V1 = [0.,;V1,;0.];
n_c3 = sum(V3.*V3.*dx);
V3 = -1./sqrt(n_c3).*V3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
figure(1)
hold on
plot(x, [0,V1',0].^2 , 'b-', 'LineWidth', 1.5 )
%
plot(x, [0,V2',0].^2 , 'g-', 'LineWidth', 1.5 )
%
plot(x, [0,V3',0].^2 , 'k-', 'LineWidth', 1.5 )
%
xlabel('$x\,(au)$','Interpreter','latex') % ,'fontsize',16
ylabel('$|\psi_{n}(x)|^{2}$','Interpreter','latex') % , 'Rotation',0
hold off
set(gca,'FontSize',16)
box on

%%%
return
end
