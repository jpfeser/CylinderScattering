function [phi_1,psi_1,phi_2,psi_2]=get_angles(wtype,cL_1,cT_1,cL_2,cT_2,phi_2,psi_2)

if strcmp(wtype,'c')
    psi_2 = asin(cT_2/cL_2*sin(phi_2));
    phi_1 = asin(cL_1/cL_2*sin(phi_2));
    psi_1 = asin(cT_1/cL_2*sin(phi_2));
elseif or(strcmp(wtype,'s_axial'),strcmp(wtype,'s_perp'))
    phi_2 = asin(cL_2/cT_2*sin(psi_2));
    phi_1 = asin(cL_1/cT_2*sin(psi_2));
    psi_1 = asin(cT_1/cT_2*sin(psi_2));
else
    fprintf('error:  to choose between +/- sign, need character assignment:\n')
    fprintf('options:  c (compression), s_axial (shear, axially aligned), s_perp (shear,perp to axis)\n')
    return
end