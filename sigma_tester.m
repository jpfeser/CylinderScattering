%% initialize
clear all
close all
PropertiesForSiGe
theta = 0;
phi = pi/2;
phi_cyl = 0;
theta_cyl = 0;
tol = 1e-12;

k = 2*pi/3e-10;
ka_vect = logspace(-3,5,100);
avect = ka_vect/k;
Psi = Get_Psi(phi,theta,phi_cyl,theta_cyl)*ones(size(ka_vect));
sigma_1 = zeros(length(avect));
sigma_2 = zeros(length(avect));
sigma_3 = zeros(length(avect));

%% case 1
for n = 1:length(avect)
    MatParams.a_NP = avect(n);
    
    p = 1;
    sigma_1(n) = GetSigmaCyl(k,p,Psi,MatParams,tol);
    
    p = 2;
    sigma_2(n) = GetSigmaCyl(k,p,Psi,MatParams,tol);
    
    p = 3;
    sigma_3(n) = GetSigmaCyl(k,p,Psi,MatParams,tol);
end

%% plots
% loglog(ka_vect,sigma_1/100e-6,'k',ka_vect,sigma_3./avect,'g',ka_vect,sigma_4./avect,'r')
loglog(ka_vect,sigma_1./(2.*avect),'b',ka_vect,sigma_2./(2.*avect),'g',ka_vect,sigma_3./(2.*avect),'r')

figure(gcf)