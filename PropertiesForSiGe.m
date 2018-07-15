
% Build Materials Properties of a Matrix / NP system
% Matrix Properties
%Silicon
MatParams.Si_vs = [8442.3 5846.2 5846.2];
MatParams.Si_omega_max = [77 28 28]*1e12;%from ioffe/Dolling(1963)
Si_A = 1.407e-19;
Si_B = 140.0177;
Si_F = 1;
Si_rho = 2329;

%Germanium
MatParams.Ge_vs = [4865.3 3566.3 3566.3];
MatParams.Ge_omega_max = [44.8 15 15]*1e12;%rad/s from ioffe/Weber (1977)
Ge_A = 2.7005e-19;
Ge_B = 64.7034;
Ge_F = 1;
Ge_rho = 5323;

%% Alloy Property Testing code % These can be used while testing for A,B & K properties for materials that make up the alloy
% MatParams.xalloy = 0.0; %alloy percentage                                                                  
% MatParams.vs = [4865.3 3566.3 3566.3]; %sound speeds                                                            
% % MatParams.A = Si_A; %Anharmonic scattering terms                                                         %change back to correct ones
% % MatParams.B = Si_B; %Anharmonic scattering terms                                                             %change back to correct ones
% MatParams.A = Ge_A;
% MatParams.B = Ge_B;
% MatParams.F = 2e-30; %alloy scattering
% MatParams.Lb = 300e-6; %boundary scattering / film thickness
% MatParams.rho = 2329; %density 

%% Alloy Code
%Alloy Properties
MatParams.xalloy = 0.5; %alloy percentage ,x: perc Ge in alloy                                                                 
MatParams.vs = MatParams.xalloy.*MatParams.Ge_vs+(1-MatParams.xalloy).*MatParams.Si_vs; %sound speeds
MatParams.omega_max = MatParams.xalloy.*MatParams.Ge_omega_max+(1-MatParams.xalloy).*MatParams.Si_omega_max; %max frequency
MatParams.kmax = MatParams.omega_max./MatParams.vs;
MatParams.A_Si = Si_A; %Anharmonic scattering terms                                                         %change back to correct ones
MatParams.B_Si = Si_B; %Anharmonic scattering terms                                                             %change back to correct ones
MatParams.A_Ge = Ge_A;
MatParams.B_Ge = Ge_B;
MatParams.F = 4.9e-30; %alloy scattering
MatParams.Lb = 300e-6; %boundary scattering / film thickness
MatParams.rho = MatParams.xalloy*Ge_rho+(1-MatParams.xalloy)*Si_rho; %density                                                                             

%% Nanoparticle Properties NiSi2
MatParams.a_NP = 10e-9; %nanocylinder radius
MatParams.vs_NP_Material = [6436.8, 3321.9, 3321.9]; %nanocylinder sound speeds
MatParams.rho_NP_Material = 4803; %nanocylinder density
MatParams.VolFrac_NP = 0.034; %volume fraction of nanocylinders                                             %change back to 0.05
MatParams.eta_NP = MatParams.VolFrac_NP/(pi*MatParams.a_NP^2); %number density (#/m2) of nanocylinders.

%% Source of the information
%%Nickel Di Silicide (NiSi2)
%Elastic properties of NiSi2, CoSi2, and PtSi by tightbinding calculations
%Giovanna Malegori and Leo Miglio PHYSICAL REVIEW B 1 OCTOBER 1993-I
%C11 = 199e9Pa 
%C44 = 53e9Pa 
%Numerical experiment on the circulation in the Japan Sea - Journal of the Oceanographical Society of Japan
%Jong-Hwan Yoon - May 1982
%rho = 4803kg/m3
%So
%c1 = 6436.8
%ci = 3321.9

%CoSi2
%Determination of the elastic constants of a cobalt disilicide
%intermetallic compound
%Guenin, Ignat and Thomas Journal of Applied Physics 68 1990
%C11 = 228e9Pa
%C44 = 83e9Pa
%Structure of Cobalt DiSilicide
%Bertaunt and Blum Seances De L Academie T231 1950 pg 627
%rho = 4940kg/m3
%So
%c1 = 6793.7
%ci = 4099

%PtSi
%Chemical bonding, elasticity, and valence force field models: A case...
%JE Klepeis, O Beckstein, O Pankratov, GLW Hart Physics Review B v64 2001
%C11 = 298.2e9Pa
%C44 = 100.1e9Pa
%Zur Kristallchemie Der B-Metall-Reichsten Phasen in Legierungen Von
%Ubergangsmetallen Der Eisentriaden Und
%K Schubert H Pfisterer Zeitschrift fur Metallkunde v41 pg 358
%rho = 12378kg/m3
%So
%c1 = 4908.3
%ci = 2843.8