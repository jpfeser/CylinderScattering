function f = pm(wtype)

if strcmp(wtype,'c')
    f = 1;
elseif or(strcmp(wtype,'s_axial'),strcmp(wtype,'s_perp'))
    f = -1;
else
    fprintf('error:  to choose between +/- sign, need character assignment:\n')
    fprintf('options:  c (compression), s_axial (shear, axially aligned), s_perp (shear,perp to axis)\n')
    return
end