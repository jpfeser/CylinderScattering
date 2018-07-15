function [Trr_inc,Trt_inc,Trz_inc,ur_inc,ut_inc,uz_inc] = incident_calc_full(n,wtype,phi,x_1,x_2,x_I,x_II,mu_1,lambda_1,mu_2,lambda_2,omega,rho_2,k_2,k_II,a)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

Cphi = cos(phi);
Sphi = sin(phi);
x_2_p = x_2*Cphi;
x_II_p = x_II*Cphi;
k2 = x_2/a;
kII= x_II/a;

C1 = lambda_2/(2*mu_2);
if n == 0
    en = 1;
else
    en = 2;
end


if strcmp(wtype,'c')
    Ka = x_2*sin(phi);
    bj = besselj(n,x_2_p);
    bjp = besselj_prime(n,x_2_p);
    Trr_inc = -2*x_2_p^2*(besselj_doubleprime(n,x_2_p)-(x_2/x_2_p)^2*C1*bj); 
    Trt_inc = 2*n*(x_2_p.*bjp- bj);
    Trz_inc = -2i*x_2_p*Ka*bjp;
    ur_inc = -x_2_p*bjp;
    ut_inc = n*bj;
    uz_inc = -1i*Ka*bj;
    %f = en*1i^n/(rho_2*omega^2);
    f = en*1i^n/(1i*k2);
elseif strcmp(wtype,'s_axial')
    Ka = x_II*sin(phi);
    bj = besselj(n,x_II_p);
    bjp = besselj_prime(n,x_II_p);
    Trr_inc = 2i*(x_II_p).^2*Sphi*besselj_doubleprime(n,x_II_p);
    Trt_inc = 2i*Sphi*n*(x_II_p.*bjp-bj);
    Trz_inc = 2*x_II_p^3/x_II*(1-0.5*(x_II/x_II_p)^2)*bjp;
    ur_inc = 1i*Ka*x_II_p/x_II*bjp;
    ut_inc = 1i*n*Ka/x_II*bj;
    uz_inc = x_II_p^2/x_II*bj;
    %f = en*1i^n/(rho_2*omega^2);
    f = en*1i^n/(1i*kII);
elseif strcmp(wtype,'s_perp')
    Ka = x_II*sin(phi);
    bj = besselj(n,x_II_p);
    bjp = besselj_prime(n,x_II_p);
    Trr_inc = -2*n*(x_II_p.*bjp-bj);
    Trt_inc = -(x_II_p).^2.*(2*besselj_doubleprime(n,x_II_p)+bj);
    Trz_inc = -1i*Ka*n*bj;
    ur_inc = -n*bj;
    ut_inc = -x_II_p.*bjp;
    uz_inc = 0;
    f = en*1i^n/(1i*kII);
    %f = en*1i^n/(rho_2*omega^2);
else
    fprintf('error:  to choose between +/- sign, need character assignment:\n')
    fprintf('options:  c (compression), s_axial (shear, axially aligned), s_perp (shear,perp to axis)\n')
    return
end

    Trr_inc = Trr_inc* f;
    Trt_inc = Trt_inc* f;
    Trz_inc = Trz_inc* f;
    ur_inc = ur_inc* f;
    ut_inc = ut_inc* f;
    uz_inc = uz_inc* f;
end

