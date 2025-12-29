%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves the one-dimensional (dim) Schrodinger equation (TISE) for the
% particle in 1D-box with length x = 0 and x = L, using finite difference method with Dirichlet boundary condition [1]. 
% Numerical solutions are compared with analytical solutions. 
%
% The atomic units are used (\hbar = m = e = 1). 
% The (TISE): H(x)*psi(x) = En*psi(x) => -0.5*d^2/dx^2*psi(x) = En*psi(x) 
% The Dirichlet boundary condition: psi(0) = psi(L) = 0.         
%
% Analytic solutions are: En_exact = pi.^2.*n.^2./(2.*L.^2), n = 1, 2, ... % exact eigenvalues
%                         psi_n(x) = sqrt(2/L)*sin(n*pi*x/L)               % exact wave function  
%
% [1] Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 29, 2025
%
%%% 2025.12.29
function [] = one_dim_box_time_independent_Schrodinger_eq
%
clc; clear one_dim_box_time_independent_Schrodinger_eq
%
% grid and initial condition
a = 0.; % x(0)
b = 5.; % x(N+1)
N = 16;  % number of grid point of x axis
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
H0_ham = -0.5.*u_mat;
%
[Vec,En] = eig(H0_ham(2:N,2:N));                                     % Eigenvalue problem with the Dirichlet boundary condition
En = diag(En);
[foo, ij] = sort(En);
En = En(ij);
[En(1),En(2),En(3),En(4),En(5)]
%
En_exact = pi.^2.*(1:5).^2./(2.*b.^2) % exact eigenvalues
%
[
%N  En(0)     En(1)     En(2)     En(3)     En(4)   
 8 0.1949    0.7498    1.5803    2.5600    3.5397
16 0.1968    0.7795    1.7258    2.9992    4.5510
32 0.1972    0.7870    1.7637    3.1179    4.8365
64 0.1974    0.7889    1.7733    3.1481    4.9101
%
%  0.1974    0.7896    1.7765    3.1583    4.9348 % <= En_exact = pi.^2.*(1:5).^2./(2.*b.^2)
];

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
V2 = -1./sqrt(n_c2).*V2;
%
V3 = Vec(:,3);                         % The unnormalized eigenfunction for the 2nd excited state,
%V1 = [0.,;V1,;0.];
n_c3 = sum(V3.*V3.*dx);
V3 = -1./sqrt(n_c3).*V3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% psi_n(x) = sqrt(2/L)*sin(n*pi*x/L) % exact wave function

psi0_exact = sqrt(2./b).*sin(1.*pi.*x./b);
psi1_exact = sqrt(2./b).*sin(2.*pi.*x./b);
psi2_exact = sqrt(2./b).*sin(3.*pi.*x./b);
%%%
figure(1)
hold on
plot(x, [0,V1',0], 'b-', 'LineWidth', 1.5 )
plot(x, psi0_exact , 'ro', 'LineWidth',1.5)
%
plot(x, [0,V2',0] + En(2), 'g-', 'LineWidth', 1.5 )
plot(x, psi1_exact + En(2), 'ro', 'LineWidth',1.5)
%
plot(x, [0,V3',0] + En(3), 'k-', 'LineWidth', 1.5 )
plot(x, psi2_exact + En(3), 'ro', 'LineWidth',1.5)
%
xlabel('$x\,(au)$','Interpreter','latex') % ,'fontsize',16
ylabel('$\psi_{n}(x)$','Interpreter','latex') % , 'Rotation',0
hold off
set(gca,'FontSize',16)
box on

%%%
return
end
