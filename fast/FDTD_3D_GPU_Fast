% FDTD 3D in three layer dispersive skin with CPML
% Author:shayan dodge
% Email address:dodgeshayan@gmail.com
%% initialize the matlab workspace
close all;
clear all;
clc;
format long
% profile on
% par=parpool(4)
%% some constants
mu_0 = 1.2566370614359173e-06;
eps_0= 8.8541878176203892e-12;
c=299792458.0;% speed of light
%% wave definition
amptidute=1;
waveforms.sinusoidal.frequency=1E9;
waveforms.gaussian.number_of_cells_per_wavelength=(2/7).*150;

T=1/waveforms.sinusoidal.frequency;
lambda=(c*T)/1;
omega=2*pi*waveforms.sinusoidal.frequency;
%% FDTD variables
number_of_cells_per_wavelength=(2/7).*150;%2.*50;
dx=lambda/number_of_cells_per_wavelength;
dy=lambda/number_of_cells_per_wavelength;
dz=lambda/number_of_cells_per_wavelength;
totalTime=500*T;
courant_factor=0.86609;
dt=1/(c*sqrt((1/dx^2)+(1/dy^2)+(1/dz^2)));
dt=courant_factor*dt;
number_of_time_steps=floor(totalTime/dt);

amptidute=1;
dtDivEps0DivDz=dt/eps_0/dz;
muSource=dtDivEps0DivDz*amptidute * 2.0*pi*omega;
%% boundary conditions
boundary.type_xn='cpml';
boundary.air_buffer_number_of_cells_xn=1; 
boundary.cpml_number_of_cells_xn=20;

boundary.type_xp = 'cpml';
boundary.air_buffer_number_of_cells_xp=1;
boundary.cpml_number_of_cells_xp=20;

boundary.type_yn = 'cpml';
boundary.air_buffer_number_of_cells_yn=1;
boundary.cpml_number_of_cells_yn=20;

boundary.type_yp = 'cpml';
boundary.air_buffer_number_of_cells_yp=1;
boundary.cpml_number_of_cells_yp=20;

boundary.type_zn = 'cpml';
boundary.air_buffer_number_of_cells_zn=ceil(1/dz);%20+40+40+40;
boundary.cpml_number_of_cells_zn=ceil(0.5/dz);

boundary.type_zp='cpml';
boundary.air_buffer_number_of_cells_zp=ceil(1/dz);
boundary.cpml_number_of_cells_zp=ceil(0.5/dz);

boundary.cpml_order = 4;
boundary.cpml_sigma_max = 1;
boundary.cpml_kappa_max = 15;
boundary.cpml_alpha_order = 1; 
boundary.cpml_alpha_max = 0.24;
boundary.cpml_eps_R= 1;
%% materialtype
%here define and initialize the arrays of material types
%air
material_type(1).eps_r=1;
material_type(1).mu_r=1;
material_type(1).sigma_e=0;
material_type(1).sigma_m=1e-20;
material_type(1).color=[1 1 1];

% ************************************
% ************************************
% Contact me to see the complete code.
% Shayan Dodge 
% dodgeshayan@gmail.com
% ************************************
% ************************************

tic
%% run_fdtd_time_marching_loop
for number=1:800

[Hx,Hy,Hz,Psi_hyz_zn,Psi_hyz_zp,Psi_hzx_xn,Psi_hyx_xp,Psi_hzy_yn,Psi_hxy_yp] = updateH(Hx,Hy,Hz,Ex,Ey,Ez,[Ex(:,1:end-1,:),AA1],[Ex(:,2:end,:),AA1],...
cat(3,Ex(:,:,1:end-1),AA2),cat(3,Ex(:,:,2:end),AA2),([Ey(1:end-1,:,:);AA4]),([Ey(2:end,:,:);AA4]),cat(3,Ey(:,:,1:end-1),AA2),cat(3,Ey(:,:,2:end),AA2)...
,[Ez(1:end-1,:,:);AA4],[Ez(2:end,:,:);AA4],([Ez(:,1:end-1,:),AA1]),([Ez(:,2:end,:),AA1]),Psi_hyx_xn,cpml_b_mx_xn,...
cpml_b1_mx_xn,cpml_a_mx_xn,cpml_a1_mx_xn,Psi_hzx_xn,CPsi_hyx_xn,CPsi_hzx_xn,...
Psi_hyx_xp,cpml_b_mx_xp,cpml_b1_mx_xp,cpml_a_mx_xp,cpml_a1_mx_xp,CPsi_hyx_xp,...
CPsi_hzx_xp,Psi_hzx_xp,Psi_hxy_yn,cpml_b_my_yn,cpml_b1_my_yn,cpml_a_my_yn,...
cpml_a1_my_yn,Psi_hzy_yn,CPsi_hxy_yn,CPsi_hzy_yn,Psi_hxy_yp,cpml_b_my_yp,...
cpml_b1_my_yp,cpml_a_my_yp,cpml_a1_my_yp,Psi_hzy_yp,CPsi_hxy_yp,CPsi_hzy_yp,...
Psi_hxz_zn,cpml_b_mz_zn,cpml_b1_mz_zn,cpml_a_mz_zn,cpml_a1_mz_zn,Psi_hyz_zn...
,CPsi_hxz_zn,CPsi_hyz_zn,Psi_hyz_zp,cpml_b_mz_zp,cpml_b1_mz_zp,cpml_a_mz_zp,...
cpml_a1_mz_zp,Psi_hxz_zp,CPsi_hxz_zp,CPsi_hyz_zp,Chxh,Chxey,Chxez,Chyh,Chyez,...
Chyex,Chzh,Chzex,Chzey,dt,omega,i,dy,dx...
,sh_1,sh_2,sh_3);  

[Ex,Ey,Ez,Psi_exz_zn,Psi_exz_zp,Psi_eyx_xn,Psi_eyx_xp,Psi_exy_yn,Psi_exy_yp] = updateE(Hx,Hy,Hz,Ex,Ey,Ez,[AA1,Hx(:,1:end-2,:),AA1],[AA1,Hx(:,2:end-1,:),AA1],cat(3,AA2,Hx(:,:,1:end-2),AA2),cat(3,AA2,Hx(:,:,2:end-1),AA2),...
[AA4;Hy(1:end-2,:,:);AA4],[AA4;Hy(2:end-1,:,:);AA4],cat(3,AA2,Hy(:,:,1:end-2),AA2),cat(3,AA2,Hy(:,:,2:end-1),AA2),[AA4;Hz(1:end-2,:,:);AA4],[AA4;Hz(2:end-1,:,:);AA4],[AA1,Hz(:,1:end-2,:),AA1],[AA1,Hz(:,2:end-1,:),AA1]...
,Cexe,Cexhz,Cexhy,Ceye,Ceyhx,Ceyhz,Ceze,Cezhx,Cezhy,...
Psi_eyx_xn,cpml_b_ex_xn,cpml_b1_ex_xn,cpml_a_ex_xn,cpml_a1_ex_xn,Psi_ezx_xn,CPsi_eyx_xn,CPsi_ezx_xn,...
Psi_eyx_xp,cpml_b_ex_xp,cpml_b1_ex_xp,cpml_a_ex_xp,cpml_a1_ex_xp,Psi_ezx_xp,CPsi_eyx_xp,CPsi_ezx_xp,...
Psi_exy_yn,cpml_b_ey_yn,cpml_b1_ey_yn,cpml_a_ey_yn,cpml_a1_ey_yn,Psi_ezy_yn,CPsi_exy_yn,CPsi_ezy_yn,...
Psi_exy_yp,cpml_b_ey_yp,cpml_b1_ey_yp,cpml_a_ey_yp,cpml_a1_ey_yp,CPsi_exy_yp,CPsi_ezy_yp,Psi_ezy_yp,...
Psi_exz_zn,cpml_b_ez_zn,cpml_b1_ez_zn,cpml_a_ez_zn,cpml_a1_ez_zn,Psi_eyz_zn,CPsi_exz_zn,CPsi_eyz_zn,...
Psi_exz_zp,cpml_b_ez_zp,cpml_b1_ez_zp,cpml_a_ez_zp,cpml_a1_ez_zp,Psi_eyz_zp,CPsi_exz_zp,CPsi_eyz_zp,...
 dx,dy,dz,SH_1,SH_2,SH_3,omega,number,dt,source,frequency,t_0,tau);  

 end
toc
% delete(par);
% profile report
