function sigma = GetSigmaCyl(k,p,Psi,MatParams, relerr_tol)
% this is for the whole thing in 2D (i.e. the wave travels in the plane
% perpendicular to the axis of the cylinder)

if nargin<5
    relerr_tol = 1e-8;
end

%unpack some things
a = MatParams.a_NP;
vs=MatParams.vs;
rho_1 = MatParams.rho_NP_Material;
vs_NP = MatParams.vs_NP_Material;
cL_1 = vs_NP(1);
cT_1 = vs_NP(2);

rho_2 = MatParams.rho;
cL_2 = vs(1);
cT_2 = vs(2);


w = warning('off','MATLAB:nearlySingularMatrix');
w = warning('off','MATLAB:illConditionedMatrix');
w = warning('off','MATLAB:SingularMatrix');

C11_1 = cL_1^2*rho_1;
C44_1 = cT_1^2*rho_1;
mu_1 = C44_1;
lambda_1 = C11_1 - 2*mu_1;

C11_2 = cL_2^2*rho_2;
C44_2 = cT_2^2*rho_2;
mu_2 = C44_2;
lambda_2 = C11_2 - 2*mu_2;

if p==1
    wtype = 'c';
    phi_2 = Psi; % for longitudinal
    psi_2 = 0; % for transverse
elseif p==2
    wtype = 's_axial';
    phi_2 = 0; % for longitudinal
    psi_2 = Psi; % for transverse
elseif p==3
    wtype = 's_perp';
    phi_2 = 0; % for longitudinal
    psi_2 = Psi; % for transverse
else
    fprintf('p is not between 1 and 3...not a valid polarization \n')
end


[phi_1_mat,psi_1_mat,phi_2_mat,psi_2_mat]=get_angles(wtype,cL_1,cT_1,cL_2,cT_2,phi_2,psi_2);
if strcmp(wtype,'c')
    phi_mat = phi_2_mat;
else
    phi_mat = psi_2_mat;
end
phi_1 = phi_1_mat(:);
psi_1 = psi_1_mat(:);
phi_2 = phi_2_mat(:);
psi_2 = psi_2_mat(:);
phi = phi_mat(:);


omega_mat = vs(p)*k;
omega_vect = omega_mat(:);

QC= zeros(1,length(omega_vect));
QS1 = QC;
QS2 = QC;
Qtot = QC;
kvect = k(:);

parfor nomega = 1:length(omega_vect)
    w = warning('off','MATLAB:nearlySingularMatrix');
w = warning('off','MATLAB:illConditionedMatrix');
w = warning('off','MATLAB:SingularMatrix');
    omega = omega_vect(nomega);
    counter_omega = nomega;
    if omega~=0
        
        %omega = 1e12;
        k_1  = omega/cL_1;
        k_2  = omega/cL_2;
        k_I  = omega/cT_1;
        k_II = omega/cT_2;
        %    a = 5e-9;
        
        x_1 = k_1*a; %internal compression wave
        x_2 = k_2*a; %external
        x_I = k_I*a;
        x_II = k_II*a;
        Ka = k_2*a*sin(phi_2(nomega));
        
        x_1_p = x_1*cos(phi_1(nomega));
        x_2_p = x_2*cos(phi_2(nomega));
        x_I_p = x_I*cos(psi_1(nomega));
        x_II_p = x_II*cos(psi_2(nomega));
        
        
        counter = 0;
        AA = zeros(1,200);
        BB = AA;
        CC = AA;
        DD = AA;
        EE = AA;
        FF = AA;
        relerr=Inf;
        Qtot_temp = 0;
        A = zeros(6,6);
        b = zeros(6,1);
        while (relerr > relerr_tol)
            n = counter;
            counter = counter +1;
            [Trr_inc,Trt_inc,Trz_inc,ur_inc,ut_inc,uz_inc] = incident_calc_full(n,wtype,phi(nomega),x_1,x_2,x_I,x_II,mu_1,lambda_1,mu_2,lambda_2,omega,rho_2,k_2,k_II,a);
            
            pm_mem = pm(wtype);
            
            bh_x2p = besselh(n,x_2_p); %reused multiple times
            bh_xIIp = besselh(n,x_II_p);
            bhp_x2p = besselh_prime(n,x_2_p);
            bhp_xIIp = besselh_prime(n,x_II_p);
            bhdp_x2p = besselh_doubleprime(n,x_2_p);
            bhdp_xIIp = besselh_doubleprime(n,x_II_p);
            
            bj_x1p = besselj(n,x_1_p);
            bj_xIp = besselj(n,x_I_p);
            bjp_x1p = besselj_prime(n,x_1_p);
            bjp_xIp = besselj_prime(n,x_I_p);
            bjdp_x1p = besselj_doubleprime(n,x_1_p);
            bjdp_xIp = besselj_doubleprime(n,x_I_p);
            
            A(1,1) = -2*x_2_p^2*(bhdp_x2p - (x_2/x_2_p)^2*(lambda_2/(2*mu_2))*bh_x2p); % An, Trr (external)
            A(1,2) = pm_mem*2*n*(x_II_p*bhp_xIIp-bh_xIIp); % Bn, Trr (external)
            A(1,3) = 2i*Ka./x_II.*(x_II_p).^2.*bhdp_xIIp; % Cn, Trr (external)
            A(1,4) = -mu_1/mu_2 * (-2)*x_1_p^2*(bjdp_x1p - (x_1/x_1_p)^2*(lambda_1/(2*mu_1))*bj_x1p); % Dn, -Trr (internal)
            A(1,5) = -mu_1/mu_2*(pm_mem*2*n*(x_I_p*bjp_xIp-bj_xIp)); %En, -Trr (internal)
            A(1,6) = -mu_1/mu_2*2i*Ka./x_I.*(x_I_p).^2.*bjdp_xIp; % Fn, -Trr (internal)
            
            A(2,1) = pm_mem*2*n*(x_2_p*bhp_x2p-bh_x2p); %An
            A(2,2) = -x_II_p^2*(2*bhdp_xIIp + bh_xIIp);
            A(2,3) = -pm_mem*2i*Ka./x_II.*n.*(x_II_p.*bhp_xIIp-bh_xIIp); %
            A(2,4) = -mu_1/mu_2*(pm_mem*2*n*(x_1_p*bjp_x1p-bj_x1p));
            A(2,5) = -mu_1/mu_2*(-x_I_p^2*(2*bjdp_xIp + bj_xIp));
            A(2,6) = -mu_1/mu_2.*(-pm_mem*2i*Ka./x_I.*n.*(x_I_p.*bjp_xIp-bj_xIp));
            
            A(3,1) =              -2i*Ka*(x_2_p).*bhp_x2p;
            A(3,2) =             pm_mem*1i*Ka.*n.*bh_xIIp;
            A(3,3) =             2*x_II_p^3/x_II*(1-0.5*(x_II/x_II_p)^2)*bhp_xIIp;
            A(3,4) = -mu_1/mu_2*(-2i*Ka*(x_1_p).*bjp_x1p);
            A(3,5) = -mu_1/mu_2*(pm_mem*1i*Ka.*n.*bj_xIp);
            A(3,6) = -mu_1/mu_2*(2*x_I_p^3/x_I*(1-0.5*(x_I/x_I_p)^2)*bjp_xIp);
            
            A(4,1) = -x_2_p.*bhp_x2p;
            A(4,2) = pm_mem*n*bh_xIIp;
            A(4,3) = 1i*Ka.*x_II_p/x_II*bhp_xIIp;
            A(4,4) =  x_1_p.*bjp_x1p;
            A(4,5) = -pm_mem*n*bj_xIp;
            A(4,6) = -1i*Ka.*x_I_p/x_I*bjp_xIp;
            
            A(5,1) = pm_mem*n*bh_x2p;
            A(5,2) = -x_II_p.*bhp_xIIp;
            A(5,3) = -pm_mem*1i*(Ka./x_II).*n.*bh_xIIp;
            A(5,4) = -pm_mem*n*bj_x1p;
            A(5,5) = x_I_p.*bjp_xIp;
            A(5,6) = pm_mem*1i*(Ka./x_I).*n.*bj_xIp;
            
            A(6,1) = -1i*Ka.*bh_x2p; % **
            A(6,2) = 0;
            A(6,3) = x_II_p^2/x_II.*bh_xIIp; %**
            A(6,4) = 1i*Ka.*bj_x1p;
            A(6,5) = 0;
            A(6,6) = -x_I_p^2/x_I.*bj_xIp;
            
            %         figure(1)
            %         spy(A)
            %        nnz(A)
            b(1,1) = -Trr_inc;
            b(2,1) = -Trt_inc;
            b(3,1) = -Trz_inc;
            b(4,1) = -ur_inc;
            b(5,1) = -ut_inc;
            b(6,1) = -uz_inc;
            
            sol = A\b;
            
            AA(counter) = sol(1);
            BB(counter) = sol(2);
            CC(counter) = sol(3);
            DD(counter) = sol(4);
            EE(counter) = sol(5);
            FF(counter) = sol(6);
            
            QC(nomega)   = 2*kvect(nomega)*(abs(AA(1))^2+(AA(:))'*AA(:));
            if imag(phi_2(nomega))~=0
                QC(nomega)=0;
            end
            QS1(nomega)  = 2*kvect(nomega)*cos(psi_2(nomega))^2*(abs(BB(1))^2+(BB(:))'*BB(:));
            QS2(nomega)  = 2*kvect(nomega)*cos(psi_2(nomega))^2*(abs(CC(1))^2+(CC(:))'*CC(:));
            Qtot(nomega) = QC(nomega) + QS1(nomega) + QS2(nomega);

            relerr = abs(Qtot(nomega) - Qtot_temp)/Qtot_temp;
            
            Qtot_temp = Qtot(nomega);
        end
        
        AA(isnan(AA))=0;
        BB(isnan(BB))=0;
        CC(isnan(CC))=0;
    else
        QC(nomega)   = 0;
        QS1(nomega)  = 0;
        QS2(nomega)  = 0;
    end
    end
    %Qtot = QC + QS1 + QS2;
    %%repackage result from vector -> matrix
    [n,m]=size(k);
    sigma = reshape(Qtot,n,m);   
end
