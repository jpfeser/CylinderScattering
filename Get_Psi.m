function Psi = Get_Psi(phi,theta,phi_cyl,theta_cyl)

hatk_x = sin(phi).*cos(theta);
hatk_y = sin(phi).*sin(theta);
hatk_z = cos(phi);

hata_x = sin(phi_cyl).*cos(theta_cyl);
hata_y = sin(phi_cyl).*sin(theta_cyl);
hata_z = cos(phi_cyl);

hatk_dot_hata = hatk_x.*hata_x + hatk_y.*hata_y +hatk_z.*hata_z;

Psi = asin(hatk_dot_hata);
end