function [Ex,Ey,Ez,Psi_exz_zn,Psi_exz_zp,Psi_eyx_xn,Psi_eyx_xp,Psi_exy_yn,Psi_exy_yp] = updateE(Hx,Hy,Hz,Ex,Ey,Ez,hx_2,hx_22,hx_3,hx_32,hy_1,hy_12,hy_3,hy_32,hz_1,hz_12,hz_2...
    ,hz_22,Cexe,Cexhz,Cexhy,Ceye,Ceyhx,Ceyhz,Ceze,Cezhx,Cezhy,...
    Psi_eyx_xn,cpml_b_ex_xn,cpml_b1_ex_xn,cpml_a_ex_xn,cpml_a1_ex_xn,Psi_ezx_xn,CPsi_eyx_xn,CPsi_ezx_xn,...
    Psi_eyx_xp,cpml_b_ex_xp,cpml_b1_ex_xp,cpml_a_ex_xp,cpml_a1_ex_xp,Psi_ezx_xp,CPsi_eyx_xp,CPsi_ezx_xp,...
    Psi_exy_yn,cpml_b_ey_yn,cpml_b1_ey_yn,cpml_a_ey_yn,cpml_a1_ey_yn,Psi_ezy_yn,CPsi_exy_yn,CPsi_ezy_yn,...
    Psi_exy_yp,cpml_b_ey_yp,cpml_b1_ey_yp,cpml_a_ey_yp,cpml_a1_ey_yp,CPsi_exy_yp,CPsi_ezy_yp,Psi_ezy_yp,...
    Psi_exz_zn,cpml_b_ez_zn,cpml_b1_ez_zn,cpml_a_ez_zn,cpml_a1_ez_zn,Psi_eyz_zn,CPsi_exz_zn,CPsi_eyz_zn,...
    Psi_exz_zp,cpml_b_ez_zp,cpml_b1_ez_zp,cpml_a_ez_zp,cpml_a1_ez_zp,Psi_eyz_zp,CPsi_exz_zp,CPsi_eyz_zp,...
    dx,dy,dz,SH_1,SH_2,SH_3,omega,number,dt,source,frequency,t_0,tau)

Psi_eyx_xn= cpml_b_ex_xn.*Psi_eyx_xn+cpml_a_ex_xn.*...
    (hz_12-hz_1);
Psi_ezx_xn= cpml_b1_ex_xn.*Psi_ezx_xn+cpml_a1_ex_xn.*...
    (hy_12-hy_1);

Ey = Ey +CPsi_eyx_xn.*Psi_eyx_xn;
Ez = Ez +CPsi_ezx_xn.*Psi_ezx_xn; 

Psi_eyx_xp=cpml_b_ex_xp.*Psi_eyx_xp+cpml_a_ex_xp.*...
    (hz_12-hz_1);
Psi_ezx_xp=cpml_b1_ex_xp.*Psi_ezx_xp+cpml_a1_ex_xp.*...
    (hy_12-hy_1);

Ey=Ey+CPsi_eyx_xp.*Psi_eyx_xp;
Ez=Ez+CPsi_ezx_xp.*Psi_ezx_xp;

Psi_exy_yn =cpml_b_ey_yn.*Psi_exy_yn+cpml_a_ey_yn.*...
    (hz_22 - hz_2);
Psi_ezy_yn=cpml_b1_ey_yn.*Psi_ezy_yn+cpml_a1_ey_yn.*...
    (hx_22-hx_2);

Ex=Ex+ CPsi_exy_yn .* Psi_exy_yn;
Ez=Ez+ CPsi_ezy_yn .* Psi_ezy_yn; 

Psi_exy_yp=cpml_b_ey_yp.*Psi_exy_yp+cpml_a_ey_yp.*...
    (hz_22-hz_2);
Psi_ezy_yp= cpml_b1_ey_yp.*Psi_ezy_yp+cpml_a1_ey_yp.*...
    (hx_22 - hx_2);

Ex=Ex+CPsi_exy_yp.*Psi_exy_yp;
Ez=Ez+CPsi_ezy_yp.*Psi_ezy_yp;

Psi_exz_zn =  cpml_b_ez_zn.*Psi_exz_zn+cpml_a_ez_zn.*...
     (hy_32- hy_3);
Psi_eyz_zn = cpml_b1_ez_zn.*Psi_eyz_zn+cpml_a1_ez_zn.*...
    (hx_32 - hx_3);

Ex  =Ex+ CPsi_exz_zn .* Psi_exz_zn;
Ey = Ey+ CPsi_eyz_zn .* Psi_eyz_zn; 

Psi_exz_zp=cpml_b_ez_zp.*Psi_exz_zp+cpml_a_ez_zp.*...
     (hy_32 - hy_3);
Psi_eyz_zp=cpml_b1_ez_zp.*Psi_eyz_zp+cpml_a1_ez_zp.*...
    (hx_32 - hx_3);

Ex=Ex+CPsi_exz_zp.*Psi_exz_zp;
Ey=Ey+CPsi_eyz_zp.*Psi_eyz_zp;
 

Ex=(Cexe.*Ex+Cexhz.*(hz_22-hz_2)+Cexhy.*(hy_32-hy_3)).*SH_1;
Ey=(Ceye.*Ey+Ceyhx.*(hx_32-hx_3)+Ceyhz.*(hz_12-hz_1)).*SH_2;
Ez=(Ceze.*Ez+Cezhy.*(hy_12-hy_1)+Cezhx.*(hx_22-hx_2)).*SH_3;

Ex=Ex+source.*exp(-((number.*dt-t_0+0.0075E-9)./tau).^2);
% sin(2*pi*frequency*number.*dt);
% exp(-((number.*dt-t_0)./tau).^2);
% -(sqrt(2*exp(1))/tau)*(number.*dt-t_0).*exp(-((number.*dt-t_0)/tau).^2);
%cos(2*pi*frequency*(number.*dt-t_0)).*exp(-((number.*dt-t_0)/tau).^2);

  
end

