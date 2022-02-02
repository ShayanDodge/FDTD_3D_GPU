% FDTD 3D in three layer non dispersive skin with CPML
% Author:shayan dodge
% Email address:dodgeshayan@gmail.com
%% initialize the matlab workspace
close all;
clear all;
clc;
%% some constants
mu_0 = 1.2566370614359173e-06;
eps_0= 8.8541878176203892e-12;
c=299792458.0;% speed of light
%% wave definition
amptidute=4.8*10^-12./1.6;
waveforms.sinusoidal.frequency=10E9;
T=1/waveforms.sinusoidal.frequency;
tau=0.4.*T;
t0=0.5*tau;
lambda=(c*T)/1.8;
omega=2*pi*waveforms.sinusoidal.frequency;
%% FDTD variables
number_of_cells_per_wavelength=180;
dx=lambda/number_of_cells_per_wavelength;
dy=lambda/number_of_cells_per_wavelength;
dz=lambda/number_of_cells_per_wavelength;
totalTime=500*T;
courant_factor=0.88;
dt=1/(c*sqrt((1/dx^2)+(1/dy^2)+(1/dz^2)));
dt=courant_factor*dt;
totalTimeStep=floor(totalTime/dt);
%% boundary conditions
boundary.type_xn='cpml';
boundary.air_buffer_number_of_cells_xn=25; 
boundary.cpml_number_of_cells_xn=-20;

boundary.type_xp = 'cpml';
boundary.air_buffer_number_of_cells_xp=25;
boundary.cpml_number_of_cells_xp=-20;

boundary.type_yn = 'cpml';
boundary.air_buffer_number_of_cells_yn=25;
boundary.cpml_number_of_cells_yn=-20;

boundary.type_yp = 'cpml';
boundary.air_buffer_number_of_cells_yp=25;
boundary.cpml_number_of_cells_yp=-20;

boundary.type_zn = 'cpml';
boundary.air_buffer_number_of_cells_zn=20;
boundary.cpml_number_of_cells_zn=-20;

boundary.type_zp='cpml';
boundary.air_buffer_number_of_cells_zp=20+40+40+40+40+40+40+40+40+40+200+200+200+200+200+200+200+200;%+200+200+200+200+200+200+200+200+200+200
boundary.cpml_number_of_cells_zp=-20;

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
%a dielectric_1
material_type(2).eps_r=2.5;
material_type(2).mu_r=1;
material_type(2).sigma_e=0;
material_type(2).sigma_m=1e-20;
material_type(2).color=[1 1 1];
%a dielectric_2
material_type(3).eps_r=5.5;
material_type(3).mu_r=1;
material_type(3).sigma_e=0;
material_type(3).sigma_m=1e-20;
material_type(3).color=[1 1 1];
%a dielectric_3
material_type(4).eps_r=59;
material_type(4).mu_r=1;
material_type(4).sigma_e=0;
material_type(4).sigma_m=1e-20;
material_type(4).color=[1 1 1];
%a dielectric_4
material_type(5).eps_r=20;
material_type(5).mu_r=1;
material_type(5).sigma_e=0;
material_type(5).sigma_m=1e-20;
material_type(5).color=[1 1 1];
%indices of material types defining air, pec, and pmc
material_type_index_air=1;
material_type_index_dielectric_1=2;
material_type_index_dielectric_2=3;
material_type_index_dielectric_3=4;
material_type_index_dielectric_4=5;
%% define_geometry
%define a brick
brick(1).min_x= 0;
brick(1).min_y= 0;
brick(1).min_z= 0;
brick(1).max_x= 0;
brick(1).max_y= 0;
brick(1).max_z= 0;
brick(1).material_type=1;

% brick(2).min_x=-10.1e-3;
% brick(2).min_y=-10.1e-3;
% brick(2).min_z= 1e-3;
% brick(2).max_x= 10.1e-3;
% brick(2).max_y= 10.1e-3;
% brick(2).max_z= 1.5e-3;
% brick(2).material_type=3;
% 
% brick(3).min_x=-10.1e-3;
% brick(3).min_y=-10.1e-3;
% brick(3).min_z= 1.5e-3;
% brick(3).max_x= 10.1e-3;
% brick(3).max_y= 10.1e-3;
% brick(3).max_z= 6.5e-3;
% brick(3).material_type=4;
% 
% brick(4).min_x=-10.1e-3;
% brick(4).min_y=-10.1e-3;
% brick(4).min_z= 6.5e-3;
% brick(4).max_x= 10.1e-3;
% brick(4).max_y= 10.1e-3;
% brick(4).max_z= 10.5e-3;
% brick(4).material_type=5;
%% calculating the domain size
number_of_brick=size(brick,2);
%find the minimum and maximum coordinates of a box encapsulating the object
number_of_objects=1;
for i=1:number_of_brick
    min_x(number_of_objects)=brick(i).min_x;
    min_y(number_of_objects)=brick(i).min_y;
    min_z(number_of_objects)=brick(i).min_z;
    max_x(number_of_objects)=brick(i).max_x;
    max_y(number_of_objects)=brick(i).max_y; 
    max_z(number_of_objects)=brick(i).max_z;
    number_of_objects=number_of_objects+1;
end
fdtd_domain.min_x=min(min_x);
fdtd_domain.min_y=min(min_y);    
fdtd_domain.min_z=min(min_z);    
fdtd_domain.max_x=max(max_x);    
fdtd_domain.max_y=max(max_y);    
fdtd_domain.max_z=max(max_z);     
% Determine the problem space boundaries including air buffers
fdtd_domain.min_x = fdtd_domain.min_x-dx *...
    boundary.air_buffer_number_of_cells_xn;
fdtd_domain.min_y = fdtd_domain.min_y-dy *...
    boundary.air_buffer_number_of_cells_yn ;
fdtd_domain.min_z = fdtd_domain.min_z-dz*...
    boundary.air_buffer_number_of_cells_zn ; 
fdtd_domain.max_x = fdtd_domain.max_x+dx *...
    boundary.air_buffer_number_of_cells_xp; 
fdtd_domain.max_y = fdtd_domain. max_y+dy *...
    boundary.air_buffer_number_of_cells_yp ;
fdtd_domain. max_z= fdtd_domain. max_z+dz *...
    boundary.air_buffer_number_of_cells_zp;
% Determine the problem space boundaries including cpml layers
if strcmp (boundary.type_xn,'cpml') &&(boundary.cpml_number_of_cells_xn>0)
    fdtd_domain.min_x = fdtd_domain.min_x -dx *...
        boundary.cpml_number_of_cells_xn;
end
if strcmp (boundary.type_xp, 'cpml') &&(boundary.cpml_number_of_cells_xp >0)
    fdtd_domain.max_x = fdtd_domain. max_x + dx *...
        boundary.cpml_number_of_cells_xp;
end
if strcmp( boundary.type_yn, 'cpml') &&(boundary.cpml_number_of_cells_yn >0)
    fdtd_domain .min_y = fdtd_domain. min_y - dy *...
        boundary.cpml_number_of_cells_yn;
end
if strcmp (boundary.type_yp, 'cpml') &&(boundary.cpml_number_of_cells_yp >0)
   fdtd_domain. max_y = fdtd_domain. max_y+ dy *...
       boundary.cpml_number_of_cells_yp ;
end
if strcmp( boundary.type_zn, 'cpml') &&(boundary.cpml_number_of_cells_zn >0)
    fdtd_domain. min_z = fdtd_domain.min_z - dz *...
        boundary.cpml_number_of_cells_zn ;
end
if strcmp (boundary.type_zp, 'cpml') && (boundary.cpml_number_of_cells_zp>0) 
    fdtd_domain. max_z = fdtd_domain. max_z + dz*...
        boundary.cpml_number_of_cells_zp;
end
%detemining the problem space size
fdtd_domain.size_x=fdtd_domain.max_x-fdtd_domain.min_x;
fdtd_domain.size_y=fdtd_domain.max_y-fdtd_domain.min_y;
fdtd_domain.size_z=fdtd_domain.max_z-fdtd_domain.min_z;
%number of cells in x, y, and z directions
nx=round(fdtd_domain.size_x/dx);
ny=round(fdtd_domain.size_y/dy);
nz=round(fdtd_domain.size_z/dz);
%adjust domain size by snapping to cells
fdtd_domain.size_x=nx*dx;
fdtd_domain.size_y=ny*dy;
fdtd_domain.size_z=nz*dz;
fdtd_domain.max_x=fdtd_domain.min_x+fdtd_domain.size_x;
fdtd_domain.max_y=fdtd_domain.min_y+fdtd_domain.size_y;
fdtd_domain.max_z=fdtd_domain.min_z+fdtd_domain.size_z;
%some frequently used auxiliary parametrs
nxp1=nx+1;  nyp1=ny+1;  nzp1=nz+1;
nxm1=nx-1;  nym1=ny-1;  nzm1=nz-1;
nxm2=nx-2;  nym2=ny-2;  nzm2=nz-2;

xcoor=linspace (fdtd_domain. min_x, fdtd_domain.max_x, nxp1);
ycoor=linspace (fdtd_domain. min_y, fdtd_domain.max_y, nyp1);
zcoor=linspace (fdtd_domain. min_z, fdtd_domain.max_z, nzp1);
%% material_3d_space
material_3d_space=ones(nx,ny,nz);
%% creating_brick
%creat the 3d object in problem space by
for ind=1:number_of_brick
%%convert brick end coordinates to node indices
blx=round((brick(ind).min_x-fdtd_domain.min_x)/dx)+1;
bly=round((brick(ind).min_y-fdtd_domain.min_y)/dy)+1;
blz=round((brick(ind).min_z-fdtd_domain.min_z)/dz)+1;
bux=round((brick(ind).max_x-fdtd_domain.min_x)/dx)+1;
buy=round((brick(ind).max_y-fdtd_domain.min_y)/dy)+1;
buz=round((brick(ind).max_z-fdtd_domain.min_z)/dz)+1;
%%assign material type of brick to the cells
material_3d_space(blx:bux-1,bly:buy-1,blz:buz-1)=brick(ind).material_type;
end
eps_r_x=ones(nx,nyp1,nzp1);
eps_r_y=ones(nxp1,ny,nzp1);
eps_r_z=ones(nxp1,nyp1,nz);
mu_r_x=ones(nxp1,ny,nz);
mu_r_y=ones(nx,nyp1,nz);
mu_r_z=ones(nx,ny,nzp1);
sigma_e_x=zeros(nx,nyp1,nzp1);
sigma_e_y=zeros(nxp1,ny,nzp1);
sigma_e_z=zeros(nxp1,nyp1,nz);
sigma_m_x=zeros(nxp1,ny,nz);
sigma_m_y=zeros(nx,nyp1,nz);
sigma_m_z=zeros(nx,ny,nzp1);
%% calculate_material_component_values
%calculate material component values by averaging
for ind=1:size(material_type,2)
    t_eps_r(ind)=material_type(ind).eps_r;
    t_mu_r(ind)=material_type(ind).mu_r;
    t_sigma_e(ind)=material_type(ind).sigma_e;
    t_sigma_m(ind)=material_type(ind).sigma_m;
end
%assign negligibly small values to t_mu_r and t_sigma_m where they are zero
%in order to prevent division by zerro error
t_mu_r(find(t_mu_r==0))=1e-20;
t_sigma_m(find(t_sigma_m==0.0000))=1e-20;
disp('calculating_eps_r_x');
%eps_r_x(i,j,k)is average of four cells(i,j,k)(i,j-1,k)(i,j,k-1)(i,j-1,k-1)
eps_r_x(1:nx,2:ny,2:nz)=0.25*(t_eps_r(material_3d_space(1:nx,2:ny,2:nz))+...
    t_eps_r(material_3d_space(1:nx,1:ny-1,2:nz))+...
    t_eps_r(material_3d_space(1:nx,2:ny,1:nz-1))+...
    t_eps_r(material_3d_space(1:nx,1:ny-1,1:nz-1)));
disp('calculating_eps_r_y');
%%eps_r_y(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j,k-1)(i-1,j,k-1)
eps_r_y(2:nx,1:ny,2:nz)=0.25*(t_eps_r(material_3d_space(2:nx,1:ny,2:nz))+...
    t_eps_r(material_3d_space(1:nx-1,1:ny,2:nz))+...
    t_eps_r(material_3d_space(2:nx,1:ny,1:nz-1))+...
    t_eps_r(material_3d_space(1:nx-1,1:ny,1:nz-1)));
disp('calculating_eps_r_z');
%eps_r_z(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j-1,k)(i-1,j-1,k)
eps_r_z(2:nx,2:ny,1:nz)=0.25*(t_eps_r(material_3d_space(2:nx,1:ny-1,1:nz))+...
    t_eps_r(material_3d_space(1:nx-1,2:ny,1:nz))+...
    t_eps_r(material_3d_space(2:nx,2:ny,1:nz))+...
    t_eps_r(material_3d_space(1:nx-1,1:ny-1,1:nz)));
disp('calculating_sigma_e_x');
%%sigma_e_x(i,j,k)is average of four cells(i,j,k)(i,j-1,k)(i,j,k-1)(i,j-1,k-1)
sigma_e_x(1:nx,2:ny,2:nz)=0.25*(t_sigma_e(material_3d_space(1:nx,2:ny,2:nz))+...
    t_sigma_e(material_3d_space(1:nx,1:ny-1,2:nz))+...
    t_sigma_e(material_3d_space(1:nx,2:ny,1:nz-1))+...
    t_sigma_e(material_3d_space(1:nx,1:ny-1,1:nz-1)));
disp('calculating_sigma_e_y');
%%sigma_e_y(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j,k-1)(i-1,j,k-1)
sigma_e_y(2:nx,1:ny,2:nz)=0.25*(t_sigma_e(material_3d_space(2:nx,1:ny,2:nz))+...
    t_sigma_e(material_3d_space(1:nx-1,1:ny,2:nz))+...
    t_sigma_e(material_3d_space(2:nx,1:ny,1:nz-1))+...
    t_sigma_e(material_3d_space(1:nx-1,1:ny,1:nz-1)));
disp('calculating_sigma_e_z');
%sigma_e_z(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j-1,k)(i-1,j-1,k)
sigma_e_z(2:nx,2:ny,1:nz)=0.25*(t_sigma_e(material_3d_space(2:nx,1:ny-1,1:nz))+...
    t_sigma_e(material_3d_space(1:nx-1,2:ny,1:nz))+...
    t_sigma_e(material_3d_space(2:nx,2:ny,1:nz))+...
    t_sigma_e(material_3d_space(1:nx-1,1:ny-1,1:nz)));
disp('calculating_sigma_m_x');
%%sigma_m_x(i,j,k)is average of two cells(i,j,k)(i-1,j,k)
sigma_m_x(2:nx,1:ny,1:nz)=0.5*(t_sigma_m(material_3d_space(2:nx,1:ny,1:nz))+...
    t_sigma_m(material_3d_space(1:nx-1,1:ny,1:nz)));
disp('calculating_sigma_m_y');
%sigma_e_y(i,j,k)is average of two cells(i,j,k)(i,j-1,k)
sigma_m_y(1:nx,2:ny,1:nz)=0.5*(t_sigma_m(material_3d_space(1:nx,2:ny,1:nz))+...
    t_sigma_m(material_3d_space(1:nx,1:ny-1,1:nz)));
disp('calculating_sigma_m_z');
%sigma_e_z(i,j,k)is average of two cells(i,j,k)(i,j,k-1)
sigma_m_z(1:nx,1:ny,2:nz)=0.5*(t_sigma_m(material_3d_space(1:nx,1:ny,2:nz))+...
    t_sigma_m(material_3d_space(1:nx,1:ny,1:nz-1)));
disp('calculating_mu_r_x');
%mu_r_x(i,j,k)is average of two cells(i,j,k)(i-1,j,k)
mu_r_x(2:nx,1:ny,1:nz)=0.5*(t_mu_r(material_3d_space(2:nx,1:ny,1:nz))+...
    t_mu_r(material_3d_space(1:nx-1,1:ny,1:nz)));
disp('calculating_mu_r_y');
%mu_r_y(i,j,k)is average of two cells(i,j,k)(i,j-1,k)
mu_r_y(1:nx,2:ny,1:nz)=0.5*(t_mu_r(material_3d_space(1:nx,2:ny,1:nz))+...
    t_mu_r(material_3d_space(1:nx,1:ny-1,1:nz)));
disp('calculating_mu_r_z');
%mu_r_z(i,j,k)is average of two cells(i,j,k)(i,j,k-1)
mu_r_z(1:nx,1:ny,2:nz)=0.5*(t_mu_r(material_3d_space(1:nx,1:ny,2:nz))+...
    t_mu_r(material_3d_space(1:nx,1:ny,1:nz-1)));

%% creating_field_array
%create and initialize field and current arrays
Hx=zeros(nxp1,ny,nz,'gpuArray');
Hy=zeros(nx,nyp1,nz,'gpuArray');
Hz=zeros(nx,ny,nzp1,'gpuArray');
Ex=zeros(nx,nyp1,nzp1,'gpuArray');
Ey=zeros(nxp1,ny,nzp1,'gpuArray');
Ez=zeros(nxp1,nyp1,nz,'gpuArray');
%% initialize_updating_coefficients
% Coeffiecients updating Ex 
Cexe=(2*eps_r_x*eps_0-dt*sigma_e_x)./(2*eps_r_x*eps_0+dt*sigma_e_x);
Cexhz=(2*dt/dy)./(2*eps_r_x*eps_0+dt*sigma_e_x);
Cexhy=-(2*dt/dz)./(2*eps_r_x*eps_0+dt*sigma_e_x);
% Coeffiecients updating Ey 
Ceye=(2*eps_r_y*eps_0-dt*sigma_e_y)./(2*eps_r_y*eps_0+dt*sigma_e_y);
Ceyhx=(2*dt/dz)./(2*eps_r_y*eps_0+dt*sigma_e_y);
Ceyhz=-(2*dt/dx)./(2*eps_r_y*eps_0+dt*sigma_e_y);
% Coeffiecients updating Ez
Ceze =(2*eps_r_z*eps_0-dt*sigma_e_z)./(2*eps_r_z.*eps_0+dt*sigma_e_z);
Cezhy=(2*dt/dx)./(2*eps_r_z*eps_0+dt*sigma_e_z);
Cezhx=-(2*dt/dy)./(2*eps_r_z.*eps_0+dt*sigma_e_z);
%general magnetic field updating coefficients
% Coeffiecients updating Hx
Chxh =(2*mu_r_x*mu_0-dt*sigma_m_x)./(2*mu_r_x*mu_0+dt*sigma_m_x);
Chxez=-(2*dt/dy)./(2*mu_r_x*mu_0+dt*sigma_m_x);
Chxey=(2*dt/dz)./(2*mu_r_x*mu_0+dt*sigma_m_x);
% Coeffiecients updating Hy
Chyh=(2*mu_r_y*mu_0-dt*sigma_m_y)./(2*mu_r_y.*mu_0+dt*sigma_m_y);
Chyex=(-2*dt/dz)./(2*mu_r_y*mu_0+dt*sigma_m_y);
Chyez=(2*dt/dx)./(2*mu_r_y.*mu_0+dt*sigma_m_y);
% Coeffiecients updating Hz
Chzh =(2*mu_r_z*mu_0-dt*sigma_m_z)./(2*mu_r_z*mu_0+dt*sigma_m_z);
Chzey=-(2*dt/dx)./(2*mu_r_z*mu_0+dt*sigma_m_z);
Chzex=(2*dt/dy)./(2*mu_r_z*mu_0+dt*sigma_m_z);

%% initialize_boundary_conditions_3d
% define logical parameters for the conditions that will be used often
n_cpml_xn = abs( boundary.cpml_number_of_cells_xn);
n_cpml_xp = abs( boundary.cpml_number_of_cells_xp);
n_cpml_yn = abs( boundary.cpml_number_of_cells_yn);
n_cpml_yp = abs( boundary.cpml_number_of_cells_yp);
n_cpml_zn= abs ( boundary.cpml_number_of_cells_zn); 
n_cpml_zp = abs( boundary.cpml_number_of_cells_zp);
% Call CPML initialization routine if any side is CPML
% Initialize CPML boundary condition
pml_order = boundary.cpml_order;% order of the polynomial distribution 
sigma_max = boundary.cpml_sigma_max;
kappa_max = boundary.cpml_kappa_max;
alpha_order = boundary.cpml_alpha_order;
alpha_max = boundary.cpml_alpha_max;
eps_R = boundary.cpml_eps_R;
% Initialize cpml for xn region
sigma_opt=sigma_max*(n_cpml_xn+1)/(sqrt(eps_R)*150*pi*dx);
rho_e=((n_cpml_xn:-1:1)-0.75)/n_cpml_xn;
rho_m=((n_cpml_xn:-1:1)-0.25)/n_cpml_xn;
sigma_ex_xn=sigma_opt*abs(rho_e).^pml_order;
sigma_mx_xn=sigma_opt*abs(rho_m).^pml_order;
kappa_ex_xn=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mx_xn=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ex_xn=alpha_max*abs(rho_e).^pml_order;
alpha_mx_xn=alpha_max*abs(rho_m).^pml_order;

cpml_b_ex_xn=exp((-dt/eps_0)...
    *((sigma_ex_xn./kappa_ex_xn)+alpha_ex_xn));
cpml_a_ex_xn=(1/dx)*(cpml_b_ex_xn-1).*sigma_ex_xn ...
    ./(kappa_ex_xn.*(sigma_ex_xn+kappa_ex_xn.*alpha_ex_xn));
cpml_b_mx_xn=exp((-dt/eps_0)...
    *((sigma_mx_xn./kappa_mx_xn)+alpha_mx_xn));
cpml_a_mx_xn=(1/dx)*(cpml_b_mx_xn-1).*sigma_mx_xn ...
    ./(kappa_mx_xn.*(sigma_mx_xn+kappa_mx_xn.*alpha_mx_xn));

Psi_eyx_xn = zeros (n_cpml_xn , ny, nzp1);
Psi_ezx_xn = zeros (n_cpml_xn , nyp1, nz);
Psi_hyx_xn = zeros (n_cpml_xn , nyp1, nz);
Psi_hzx_xn = zeros (n_cpml_xn , ny, nzp1 );

CPsi_eyx_xn = Ceyhz (2:n_cpml_xn +1,:,:)* dx;
CPsi_ezx_xn = Cezhy (2:n_cpml_xn +1,:,:)* dx;
CPsi_hyx_xn = Chyez (1:n_cpml_xn , : , :)*dx;
CPsi_hzx_xn = Chzey (1:n_cpml_xn,:,:) * dx ;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for i = 1: n_cpml_xn
Ceyhz (i+1,:,:) = Ceyhz (i +1,:,:)/ kappa_ex_xn (i); 
Cezhy (i+1,:,:) = Cezhy(i+1,:,:)/kappa_ex_xn (i);
Chyez(i,:,:) = Chyez(i,:,:)/ kappa_mx_xn(i);
Chzey(i,:,:) = Chzey(i, :,:)/ kappa_mx_xn(i);
end
% Initialize cpml for xp region
sigma_opt=sigma_max*(n_cpml_xp+1)/(sqrt(eps_R)*150*pi*dx);
rho_e=((1:1:n_cpml_xp)-0.75)/n_cpml_xp;
rho_m=((1:1:n_cpml_xp)-0.25)/n_cpml_xp;
sigma_ex_xp=sigma_opt*abs(rho_e).^pml_order;
sigma_mx_xp=sigma_opt*abs(rho_m).^pml_order;
kappa_ex_xp=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mx_xp=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ex_xp=alpha_max*abs(rho_e).^pml_order;
alpha_mx_xp=alpha_max*abs(rho_m).^pml_order;

cpml_b_ex_xp=exp((-dt/eps_0)...
    *((sigma_ex_xp./kappa_ex_xp)+alpha_ex_xp));
cpml_a_ex_xp=(1/dx)*(cpml_b_ex_xp-1).*sigma_ex_xp ...
    ./(kappa_ex_xp.*(sigma_ex_xp+kappa_ex_xp.*alpha_ex_xp));
cpml_b_mx_xp=exp((-dt/eps_0)...
    *((sigma_mx_xp./kappa_mx_xp)+alpha_mx_xp));
cpml_a_mx_xp=1/dx*(cpml_b_mx_xp-1).*sigma_mx_xp ...
    ./(kappa_mx_xp.*(sigma_mx_xp+kappa_mx_xp.*alpha_mx_xp));

Psi_eyx_xp = zeros (n_cpml_xp , ny, nzp1);
Psi_ezx_xp = zeros (n_cpml_xp , nyp1, nz);
Psi_hyx_xp = zeros (n_cpml_xp , nyp1, nz);
Psi_hzx_xp = zeros (n_cpml_xp , ny, nzp1 );

CPsi_eyx_xp = Ceyhz (nxp1-n_cpml_xp:nx,:,:)* dx;
CPsi_ezx_xp = Cezhy (nxp1-n_cpml_xp:nx,:,:)* dx;
CPsi_hyx_xp = Chyez (nxp1-n_cpml_xp:nx,:,:)* dx;
CPsi_hzx_xp = Chzey (nxp1-n_cpml_xp:nx,:,:)* dx ;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for i = 1:n_cpml_xp
Ceyhz (nx-n_cpml_xp+i,:,:) = Ceyhz(nx-n_cpml_xp+i,:,:)/ kappa_ex_xp (i);
Cezhy (nx-n_cpml_xp+i,:,:) = Cezhy(nx-n_cpml_xp+i,:,:)/ kappa_ex_xp (i);
Chyez (nx-n_cpml_xp+i,:,:) = Chyez(nx-n_cpml_xp+i,:,:)/ kappa_mx_xp (i);
Chzey (nx-n_cpml_xp+i,:,:) = Chzey(nx-n_cpml_xp+i,:,:)/ kappa_mx_xp (i);
end

% Initialize cpml for yn region
sigma_opt=sigma_max*(n_cpml_yn+1)/(sqrt(eps_R)*150*pi*dy);
rho_e=((n_cpml_yn:-1:1)-0.75)/n_cpml_yn;
rho_m=((n_cpml_yn:-1:1)-0.25)/n_cpml_yn;
sigma_ey_yn=sigma_opt*abs(rho_e).^pml_order;
sigma_my_yn=sigma_opt*abs(rho_m).^pml_order;
kappa_ey_yn=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_my_yn=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ey_yn=alpha_max*abs(rho_e).^pml_order;
alpha_my_yn=alpha_max*abs(rho_m).^pml_order;

cpml_b_ey_yn=exp((-dt/eps_0)...
    *((sigma_ey_yn./kappa_ey_yn)+alpha_ey_yn));
cpml_a_ey_yn=1/dy*(cpml_b_ey_yn-1).*sigma_ey_yn ...
    ./(kappa_ey_yn.*(sigma_ey_yn+kappa_ey_yn.*alpha_ey_yn));
cpml_b_my_yn=exp((-dt/eps_0)...
    *((sigma_my_yn./kappa_my_yn)+alpha_my_yn));
cpml_a_my_yn=1/dy*(cpml_b_my_yn-1).*sigma_my_yn ...
    ./(kappa_my_yn.*(sigma_my_yn+kappa_my_yn.*alpha_my_yn));

Psi_exy_yn = zeros (nx,n_cpml_yn,nzp1); 
Psi_ezy_yn = zeros (nxp1,n_cpml_yn,nz);
Psi_hxy_yn = zeros (nxp1,n_cpml_yn,nz);
Psi_hzy_yn = zeros(nx,n_cpml_yn, nzp1 );
% Create and initialize 2D cpml convolution coefficients 
% Notice that Ey (1,:,:) and Ez (1,:,:) are not updated by cmpl 
CPsi_exy_yn = Cexhz (:,2:n_cpml_yn+1,:)*dy;
CPsi_ezy_yn = Cezhx (:,2:n_cpml_yn+1,:)*dy;
CPsi_hxy_yn = Chxez (:,1:n_cpml_yn  ,:)*dy;
CPsi_hzy_yn = Chzex (:,1:n_cpml_yn  ,:)*dy;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for j = 1: n_cpml_yn
Cexhz (:,j+1,:) = Cexhz (:,j+1,:)/ kappa_ey_yn (j);
Cezhx (:,j+1,:) = Cezhx (:,j+1,:)/ kappa_ey_yn (j);
Chxez (:,j  ,:) = Chxez (:,j  ,:)/ kappa_my_yn (j);
Chzex (:,j  ,:) = Chzex (:,j  ,:)/ kappa_my_yn (j);
end
% Initialize cpml for yp region
sigma_opt=sigma_max*(n_cpml_yp+1)/(sqrt(eps_R)*150*pi*dy);
rho_e=((1:1:n_cpml_yp)-0.75)/n_cpml_yp;
rho_m=((1:1:n_cpml_yp)-0.25)/n_cpml_yp;
sigma_ey_yp=sigma_opt*abs(rho_e).^pml_order;
sigma_my_yp=sigma_opt*abs(rho_m).^pml_order;
kappa_ey_yp=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_my_yp=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ey_yp=alpha_max*abs(rho_e).^pml_order;
alpha_my_yp=alpha_max*abs(rho_m).^pml_order;

cpml_b_ey_yp=exp((-dt/eps_0)...
    *((sigma_ey_yp./kappa_ey_yp)+alpha_ey_yp));
cpml_a_ey_yp=1/dx*(cpml_b_ey_yp-1).*sigma_ey_yp ...
    ./(kappa_ey_yp.*(sigma_ey_yp+kappa_ey_yp.*alpha_ey_yp));
cpml_b_my_yp=exp((-dt/eps_0)...
    *((sigma_my_yp./kappa_my_yp)+alpha_my_yp));
cpml_a_my_yp=1/dy*(cpml_b_my_yp-1).*sigma_my_yp ...
    ./(kappa_my_yp.*(sigma_my_yp+kappa_my_yp.*alpha_my_yp));

Psi_exy_yp = zeros (nx,n_cpml_yp,nzp1); 
Psi_ezy_yp = zeros (nxp1,n_cpml_yp,nz);
Psi_hxy_yp = zeros (nxp1,n_cpml_yp,nz);
Psi_hzy_yp = zeros (nx,n_cpml_yp,nzp1);
% Create and initialize 2D cpml convolution coefficients
% Notice that Ey (nxp1,:,:) and Ez (nxp1 ,:,:) are not updated by cmpl 
CPsi_exy_yp = Cexhz (:,nyp1-n_cpml_yp:ny,:)* dy;
CPsi_ezy_yp = Cezhx (:,nyp1-n_cpml_yp:ny,:)* dy;
CPsi_hxy_yp = Chxez (:,nyp1-n_cpml_yp:ny,:)* dy;
CPsi_hzy_yp = Chzex (:,nyp1-n_cpml_yp:ny,:)* dy;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (nxp1,:,:) and Ez (nxp1 ,:,:) are not updated by cmpl 
for j = 1:n_cpml_yp
Cexhz (:,ny-n_cpml_yp+j,:) = Cexhz (:,ny-n_cpml_yp+j,:)/ kappa_ey_yp (j);
Cezhx (:,ny-n_cpml_yp+j,:) = Cezhx (:,ny-n_cpml_yp+j,:)/ kappa_ey_yp (j);
Chxez (:,ny-n_cpml_yp+j,:) = Chxez (:,ny-n_cpml_yp+j,:)/ kappa_my_yp (j);
Chzex (:,ny-n_cpml_yp+j,:) = Chzex (:,ny-n_cpml_yp+j,:)/ kappa_my_yp (j);
end

% Initialize cpml for zn region
sigma_opt=sigma_max*(n_cpml_zn+1)/(sqrt(eps_R)*150*pi*dz);
rho_e=((n_cpml_zn:-1:1)-0.75)/n_cpml_zn;
rho_m=((n_cpml_zn:-1:1)-0.25)/n_cpml_zn;
sigma_ez_zn=sigma_opt*abs(rho_e).^pml_order;
sigma_mz_zn=sigma_opt*abs(rho_m).^pml_order;
kappa_ez_zn=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mz_zn=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ez_zn=alpha_max*abs(rho_e).^pml_order;
alpha_mz_zn=alpha_max*abs(rho_m).^pml_order;

cpml_b_ez_zn=exp((-dt/eps_0)...
    *((sigma_ez_zn./kappa_ez_zn)+alpha_ez_zn));
cpml_a_ez_zn=1/dz*(cpml_b_ez_zn-1).*sigma_ez_zn ...
    ./(kappa_ez_zn.*(sigma_ez_zn+kappa_ez_zn.*alpha_ez_zn));
cpml_b_mz_zn=exp((-dt/eps_0)...
    *((sigma_mz_zn./kappa_mz_zn)+alpha_mz_zn));
cpml_a_mz_zn=1/dz*(cpml_b_mz_zn-1).*sigma_mz_zn ...
    ./(kappa_mz_zn.*(sigma_mz_zn+kappa_mz_zn.*alpha_mz_zn));

Psi_eyz_zn = zeros (nxp1,ny,n_cpml_zn);
Psi_exz_zn = zeros (nx,nyp1,n_cpml_zn);
Psi_hyz_zn = zeros (nx,nyp1,n_cpml_zn);
Psi_hxz_zn = zeros (nxp1,ny,n_cpml_zn );
% Create and initialize 2D cpml convolution coefficients 
% Notice that Ey (1,:,:) and Ez (1,:,:) are not updated by cmpl 
CPsi_eyz_zn = Ceyhx (:,:,2: n_cpml_zn+1)* dz;
CPsi_exz_zn = Cexhy (:,:,2: n_cpml_zn+1)* dz;
CPsi_hyz_zn = Chyex (:,:,1: n_cpml_zn)* dz;
CPsi_hxz_zn = Chxey (:,:,1: n_cpml_zn)* dz ;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for j = 1: n_cpml_zn
Cexhy (:,:,j+1) = Cexhy (:,:,j+1)/ kappa_ez_zn (j); 
Ceyhx (:,:,j+1) = Ceyhx(:,:,j+1)/kappa_ez_zn (j);
Chxey (:,:,j) = Chxey(:,:,j)/ kappa_mz_zn(j);
Chyex (:,:,j) = Chyex(:,:,j)/ kappa_mz_zn(j);
end
% Initialize cpml for zp region
sigma_opt=sigma_max*(n_cpml_zp+1)/(sqrt(eps_R)*150*pi*dz);
rho_e=((1:1:n_cpml_zp)-0.75)/n_cpml_zp;
rho_m=((1:1:n_cpml_zp)-0.25)/n_cpml_zp;
sigma_ez_zp=sigma_opt*abs(rho_e).^pml_order;
sigma_mz_zp=sigma_opt*abs(rho_m).^pml_order;
kappa_ez_zp=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mz_zp=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ez_zp=alpha_max*abs(rho_e).^pml_order;
alpha_mz_zp=alpha_max*abs(rho_m).^pml_order;

cpml_b_ez_zp=exp((-dt/eps_0)...
    *((sigma_ez_zp./kappa_ez_zp)+alpha_ez_zp));
cpml_a_ez_zp=1/dz*(cpml_b_ez_zp-1).*sigma_ez_zp ...
    ./(kappa_ez_zp.*(sigma_ez_zp+kappa_ez_zp.*alpha_ez_zp));
cpml_b_mz_zp=exp((-dt/eps_0)...
    *((sigma_mz_zp./kappa_mz_zp)+alpha_mz_zp));
cpml_a_mz_zp=1/dz*(cpml_b_mz_zp-1).*sigma_mz_zp ...
    ./(kappa_mz_zp.*(sigma_mz_zp+kappa_mz_zp.*alpha_mz_zp));

Psi_eyz_zp = zeros (nxp1,ny,n_cpml_zp);
Psi_exz_zp = zeros (nx,nyp1,n_cpml_zp);
Psi_hyz_zp = zeros (nx,nyp1,n_cpml_zp);
Psi_hxz_zp = zeros (nxp1,ny,n_cpml_zp);
% Create and initialize 2D cpml convolution coefficients 
% Notice that Ey (1,:,:) and Ez (1,:,:) are not updated by cmpl 
CPsi_eyz_zp = Ceyhx (:,:,nzp1-n_cpml_zp:nz)* dz;
CPsi_exz_zp = Cexhy (:,:,nzp1-n_cpml_zp:nz)* dz;
CPsi_hyz_zp = Chyex (:,:,nzp1-n_cpml_zp:nz)* dz;
CPsi_hxz_zp = Chxey (:,:,nzp1-n_cpml_zp:nz)* dz;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for i = 1: n_cpml_zp
Cexhy(:,:,nz-n_cpml_zp+i) = Cexhy(:,:,nz-n_cpml_zp+i)/kappa_ez_zp(i); 
Ceyhx(:,:,nz-n_cpml_zp+i) = Ceyhx(:,:,nz-n_cpml_zp+i)/kappa_ez_zp(i);
Chxey(:,:,nz-n_cpml_zp+i) = Chxey(:,:,nz-n_cpml_zp+i)/kappa_mz_zp(i);
Chyex(:,:,nz-n_cpml_zp+i) = Chyex(:,:,nz-n_cpml_zp+i)/kappa_mz_zp(i);
end
%% ====== Gaussian Source ===========
dtDivEps0DivDz=dt/eps_0/dz;
muSource=dtDivEps0DivDz*amptidute * -2.0 /T/T;
%% initialize_output_parameters_3d
Hx_sample=zeros(nxp1-21,ny-21,nz-81)  ;
Hy_sample=zeros(nx-21,nyp1-21,nz-81)  ;
Hz_sample=zeros(nx-21,ny-21,nzp1-81)  ;
Ex_sample=zeros(nx-21,nyp1-21,nzp1-81);
Ey_sample=zeros(nxp1-21,ny-21,nzp1-81);
Ez_sample=zeros(nxp1-21,nyp1-21,nz-81);
% Hh_1=figure
% for ind=1:number_of_brick
% hold on
% patch([brick(ind).min_x brick(ind).min_x brick(ind).max_x ...
%     brick(ind).max_x],[brick(ind).min_y brick(ind).max_y...
%     brick(ind).max_y brick(ind).min_y],[brick(ind).min_z...
%     brick(ind).min_z brick(ind).min_z brick(ind).min_z],...
%     material_type(ind+1).color)
% patch([brick(ind).min_x brick(ind).min_x brick(ind).max_x...
%     brick(ind).max_x], [brick(ind).min_y brick(ind).max_y...
%     brick(ind).max_y brick(ind).min_y], [brick(ind).max_z...
%     brick(ind).max_z brick(ind).max_z brick(ind).max_z],...
%     material_type(ind+1).color)
% patch([brick(ind).min_x brick(ind).min_x brick(ind).min_x...
%     brick(ind).min_x], [brick(ind).min_y brick(ind).max_y...
%     brick(ind).max_y brick(ind).min_y], [brick(ind).min_z...
%     brick(ind).min_z brick(ind).max_z brick(ind).max_z],...
%     material_type(ind+1).color)
% patch([brick(ind).max_x brick(ind).max_x brick(ind).max_x...
%     brick(ind).max_x], [brick(ind).min_y brick(ind).max_y...
%     brick(ind).max_y brick(ind).min_y], [brick(ind).min_z...
%     brick(ind).min_z brick(ind).max_z brick(ind).max_z],...
%     material_type(ind+1).color)
% patch([brick(ind).min_x brick(ind).min_x brick(ind).max_x...
%     brick(ind).max_x], [brick(ind).min_y brick(ind).min_y...
%     brick(ind).min_y brick(ind).min_y], [brick(ind).min_z...
%     brick(ind).max_z brick(ind).max_z brick(ind).min_z],...
%     material_type(ind+1).color)
% patch([brick(ind).min_x brick(ind).min_x brick(ind).max_x...
%     brick(ind).max_x], [brick(ind).max_y brick(ind).max_y...
%     brick(ind).max_y brick(ind).max_y], [brick(ind).min_z...
%     brick(ind).max_z brick(ind).max_z brick(ind).min_z],...
%     material_type(ind+1).color)
% alpha(.009)
% view(3)
% end



% cpml_b_mx_xn=cpml_b_mx_xn';
% cpml_a_mx_xn=cpml_a_mx_xn';
% cpml_b_mx_xp=cpml_b_mx_xp';
% cpml_a_mx_xp=cpml_a_mx_xp';
% cpml_b_ex_xn=cpml_b_ex_xn';
% cpml_a_ex_xn=cpml_a_ex_xn';
% cpml_b_ex_xp=cpml_b_ex_xp';
% cpml_a_ex_xp=cpml_a_ex_xp';
% 
% cpml_b_mz_zn=permute(cpml_b_mz_zn,[1 3 2]);
% cpml_a_mz_zn=permute(cpml_a_mz_zn,[1 3 2]);
% cpml_b_mz_zp=permute(cpml_b_mz_zp,[1 3 2]);
% cpml_a_mz_zp=permute(cpml_a_mz_zp,[1 3 2]);
% 
% cpml_b_ez_zn=permute(cpml_b_ez_zn,[1 3 2]);
% cpml_a_ez_zn=permute(cpml_a_ez_zn,[1 3 2]);
% cpml_b_ez_zp=permute(cpml_b_ez_zp,[1 3 2]);
% cpml_a_ez_zp=permute(cpml_a_ez_zp,[1 3 2]);

    cpml_b_ex_xn=gpuArray(cpml_b_ex_xn);
    cpml_b_ex_xp=gpuArray(cpml_b_ex_xp);
    cpml_b_ey_yn=gpuArray(cpml_b_ey_yn);
    cpml_b_ey_yp=gpuArray(cpml_b_ey_yp);
    cpml_b_ez_zn=gpuArray(cpml_b_ez_zn);
    cpml_b_ez_zp=gpuArray(cpml_b_ez_zp);
    cpml_a_ex_xn=gpuArray(cpml_a_ex_xn);
    cpml_a_ex_xp=gpuArray(cpml_a_ex_xp);
    cpml_a_ey_yn=gpuArray(cpml_a_ey_yn);
    cpml_a_ey_yp=gpuArray(cpml_a_ey_yp);
    cpml_a_ez_zn=gpuArray(cpml_a_ez_zn);
    cpml_a_ez_zp=gpuArray(cpml_a_ez_zp);
    cpml_b_mx_xn=gpuArray(cpml_b_mx_xn);
    cpml_b_mx_xp=gpuArray(cpml_b_mx_xp);
    cpml_b_my_yn=gpuArray(cpml_b_my_yn);
    cpml_b_my_yp=gpuArray(cpml_b_my_yp);
    cpml_b_mz_zn=gpuArray(cpml_b_mz_zn);
    cpml_b_mz_zp=gpuArray(cpml_b_mz_zp);
    cpml_a_mx_xn=gpuArray(cpml_a_mx_xn);
    cpml_a_mx_xp=gpuArray(cpml_a_mx_xp);
    cpml_a_my_yn=gpuArray(cpml_a_my_yn);
    cpml_a_my_yp=gpuArray(cpml_a_my_yp);
    cpml_a_mz_zn=gpuArray(cpml_a_mz_zn);
    cpml_a_mz_zp=gpuArray(cpml_a_mz_zp);
%     cpml_b1_ex_xn=gpuArray(cpml_b1_ex_xn);
%     cpml_b1_ex_xp=gpuArray(cpml_b1_ex_xp);
%     cpml_b1_ey_yn=gpuArray(cpml_b1_ey_yn);
%     cpml_b1_ey_yp=gpuArray(cpml_b1_ey_yp);
%     cpml_b1_ez_zn=gpuArray(cpml_b1_ez_zn);
%     cpml_b1_ez_zp=gpuArray(cpml_b1_ez_zp);
%     cpml_a1_ex_xn=gpuArray(cpml_a1_ex_xn);
%     cpml_a1_ex_xp=gpuArray(cpml_a1_ex_xp);
%     cpml_a1_ey_yn=gpuArray(cpml_a1_ey_yn);
%     cpml_a1_ey_yp=gpuArray(cpml_a1_ey_yp);
%     cpml_a1_ez_zn=gpuArray(cpml_a1_ez_zn);
%     cpml_a1_ez_zp=gpuArray(cpml_a1_ez_zp);
%     cpml_b1_mx_xn=gpuArray(cpml_b1_mx_xn);
%     cpml_b1_mx_xp=gpuArray(cpml_b1_mx_xp);
%     cpml_b1_my_yn=gpuArray(cpml_b1_my_yn);
%     cpml_b1_my_yp=gpuArray(cpml_b1_my_yp);
%     cpml_b1_mz_zn=gpuArray(cpml_b1_mz_zn);
%     cpml_b1_mz_zp=gpuArray(cpml_b1_mz_zp);
%     cpml_a1_mx_xn=gpuArray(cpml_a1_mx_xn);
%     cpml_a1_mx_xp=gpuArray(cpml_a1_mx_xp);
%     cpml_a1_my_yn=gpuArray(cpml_a1_my_yn);
%     cpml_a1_my_yp=gpuArray(cpml_a1_my_yp);
%     cpml_a1_mz_zn=gpuArray(cpml_a1_mz_zn);
%     cpml_a1_mz_zp=gpuArray(cpml_a1_mz_zp);  

      CPsi_eyx_xn=gpuArray(CPsi_eyx_xn); 
      CPsi_eyx_xp=gpuArray(CPsi_eyx_xp); 
      CPsi_ezx_xn=gpuArray(CPsi_ezx_xn); 
      CPsi_ezx_xp=gpuArray(CPsi_ezx_xp); 
      CPsi_ezy_yn=gpuArray(CPsi_ezy_yn); 
      CPsi_ezy_yp=gpuArray(CPsi_ezy_yp); 
      CPsi_exy_yn=gpuArray(CPsi_exy_yn); 
      CPsi_exy_yp=gpuArray(CPsi_exy_yp); 
      CPsi_eyz_zn=gpuArray(CPsi_eyz_zn); 
      CPsi_eyz_zp=gpuArray(CPsi_eyz_zp); 
      CPsi_exz_zn=gpuArray(CPsi_exz_zn); 
      CPsi_exz_zp=gpuArray(CPsi_exz_zp); 
      
      CPsi_hyx_xn=gpuArray(CPsi_hyx_xn); 
      CPsi_hyx_xp=gpuArray(CPsi_hyx_xp); 
      CPsi_hzx_xn=gpuArray(CPsi_hzx_xn); 
      CPsi_hzx_xp=gpuArray(CPsi_hzx_xp); 
      CPsi_hzy_yn=gpuArray(CPsi_hzy_yn); 
      CPsi_hzy_yp=gpuArray(CPsi_hzy_yp); 
      CPsi_hxy_yn=gpuArray(CPsi_hxy_yn); 
      CPsi_hxy_yp=gpuArray(CPsi_hxy_yp); 
      CPsi_hyz_zn=gpuArray(CPsi_hyz_zn); 
      CPsi_hyz_zp=gpuArray(CPsi_hyz_zp); 
      CPsi_hxz_zn=gpuArray(CPsi_hxz_zn); 
      CPsi_hxz_zp=gpuArray(CPsi_hxz_zp);     
      
      Psi_eyx_xn=gpuArray(Psi_eyx_xn); 
      Psi_eyx_xp=gpuArray(Psi_eyx_xp); 
      Psi_ezx_xn=gpuArray(Psi_ezx_xn); 
      Psi_ezx_xp=gpuArray(Psi_ezx_xp); 
      Psi_ezy_yn=gpuArray(Psi_ezy_yn); 
      Psi_ezy_yp=gpuArray(Psi_ezy_yp); 
      Psi_exy_yn=gpuArray(Psi_exy_yn); 
      Psi_exy_yp=gpuArray(Psi_exy_yp); 
      Psi_eyz_zn=gpuArray(Psi_eyz_zn); 
      Psi_eyz_zp=gpuArray(Psi_eyz_zp); 
      Psi_exz_zn=gpuArray(Psi_exz_zn); 
      Psi_exz_zp=gpuArray(Psi_exz_zp); 
      
      Psi_hyx_xn=gpuArray(Psi_hyx_xn); 
      Psi_hyx_xp=gpuArray(Psi_hyx_xp); 
      Psi_hzx_xn=gpuArray(Psi_hzx_xn); 
      Psi_hzx_xp=gpuArray(Psi_hzx_xp); 
      Psi_hzy_yn=gpuArray(Psi_hzy_yn); 
      Psi_hzy_yp=gpuArray(Psi_hzy_yp); 
      Psi_hxy_yn=gpuArray(Psi_hxy_yn); 
      Psi_hxy_yp=gpuArray(Psi_hxy_yp); 
      Psi_hyz_zn=gpuArray(Psi_hyz_zn); 
      Psi_hyz_zp=gpuArray(Psi_hyz_zp); 
      Psi_hxz_zn=gpuArray(Psi_hxz_zn); 
      Psi_hxz_zp=gpuArray(Psi_hxz_zp);       
      
      Cexe=gpuArray(Cexe);
      Cexhz=gpuArray(Cexhz);
      Cexhy=gpuArray(Cexhy);

      Ceye=gpuArray(Ceye);
      Ceyhx=gpuArray(Ceyhx);
      Ceyhz=gpuArray(Ceyhz);

     Ceze=gpuArray(Ceze);
     Cezhy=gpuArray(Cezhy);
     Cezhx=gpuArray(Cezhx);

     Chxh=gpuArray(Chxh);
     Chxh=gpuArray(Chxh);
     Chxey=gpuArray(Chxey);
     Chxey=gpuArray(Chxey);
     Chxez=gpuArray(Chxez);
     Chxez=gpuArray(Chxez);

     Chyh=gpuArray(Chyh);
     Chyh=gpuArray(Chyh);
     Chyex=gpuArray(Chyex);
     Chyex=gpuArray(Chyex);
     Chyez=gpuArray(Chyez);
     Chyez=gpuArray(Chyez);

     Chzh=gpuArray(Chzh);
     Chzh=gpuArray(Chzh);
     Chzex=gpuArray(Chzex);
     Chzex=gpuArray(Chzex);
     Chzey=gpuArray(Chzey);
     Chzey=gpuArray(Chzey);
tic
%% run_fdtd_time_marching_loop
for number=1:10
%     E_0P=Ex((nx+1)/2,(ny+1)/2,100);
%% update_magnetic_fields_CPML
% Psi_hyx_xn(:,:,:) = cpml_b_mx_xn.*Psi_hyx_xn+cpml_a_mx_xn.*...
%     ( Ez (2:n_cpml_xn+1,:,:) - Ez(1:n_cpml_xn,:,:));
% Psi_hzx_xn(:,:,:) = cpml_b_mx_xn.* Psi_hzx_xn+cpml_a_mx_xn.*...
%     ( Ey (2:n_cpml_xn+1,:,:) - Ey(1:n_cpml_xn ,:,:));
% Psi_hyx_xp (:,:,:) = cpml_b_mx_xp.* Psi_hyx_xp+...
%     cpml_a_mx_xp.*(Ez (nx-n_cpml_xp+2:nx+1,:,:)-Ez(nx-n_cpml_xp+1:nx,:,:));
% Psi_hzx_xp (:,:,:) = cpml_b_mx_xp.* Psi_hzx_xp+...
%     cpml_a_mx_xp.*(Ey (nx-n_cpml_xp+2:nx+1,:,:)-Ey(1+nx-n_cpml_xp:nx,:,:));
% 
% Psi_hxy_yn(:,:,:)=cpml_b_my_yn.*Psi_hxy_yn+cpml_a_my_yn.*...
%     (Ez(:,2:n_cpml_yn+1,:) - Ez(:,1:n_cpml_yn,:));
% Psi_hzy_yn(:,:,:)=cpml_b_my_yn.*Psi_hzy_yn+cpml_a_my_yn.*...
%     (Ex(:,2:n_cpml_yn+1,:) - Ex(:,1:n_cpml_yn,:));
% Psi_hxy_yp(:,:,:)=cpml_b_my_yp.*Psi_hxy_yp(:,:,:)+cpml_a_my_yp.*...
%     (Ez(:,ny-n_cpml_yp+2:ny+1,:)-Ez(:,ny-n_cpml_yp+1:ny,:));
% Psi_hzy_yp(:,:,:)=cpml_b_my_yp.*Psi_hzy_yp(:,:,:)+cpml_a_my_yp.*...
%     (Ex(:,ny-n_cpml_yp+2:ny+1,:)-Ex(:,1+ny-n_cpml_yp:ny,:));
% 
% Psi_hxz_zn(:,:,:)=(cpml_b_mz_zn).*Psi_hxz_zn(:,:,:)+(cpml_a_mz_zn).*...
%     (Ey(:,:,2:n_cpml_zn+1) - Ey(:,:,1:n_cpml_zn));
% Psi_hyz_zn(:,:,:)=(cpml_b_mz_zn).*Psi_hyz_zn(:,:,:)+(cpml_a_mz_zn).*...
%     (Ex(:,:,2:n_cpml_zn+1) - Ex(:,:,1:n_cpml_zn));
% Psi_hyz_zp (:,:,:) = (cpml_b_mz_zp).*Psi_hyz_zp(:,:,:)+(cpml_a_mz_zp).*...
%     (Ex(:,:,nz-n_cpml_zp+2:nz+1)-Ex(:,:,nz-n_cpml_zp+1:nz));
% Psi_hxz_zp (:,:,:) = (cpml_b_mz_zp).*Psi_hxz_zp(:,:,:)+(cpml_a_mz_zp).*...
%     (Ey(:,:,nz-n_cpml_zp+2:nz+1)-Ey(:,:,nz-n_cpml_zp+1:nz));
% 
% Hy(1:n_cpml_xn,:,:)=Hy(1:n_cpml_xn,:,:)+CPsi_hyx_xn(:,:,:).*Psi_hyx_xn(:,:,:);
% Hz(1:n_cpml_xn,:,:)=Hz(1:n_cpml_xn,:,:)+CPsi_hzx_xn(:,:,:).*Psi_hzx_xn(:,:,:);
% 
% Hy(nx-n_cpml_xp+1:nx,:,:)=Hy(nx-n_cpml_xp+1:nx,:,:)+CPsi_hyx_xp(:,:,:).*Psi_hyx_xp(:,:,:);
% Hz(nx-n_cpml_xp+1:nx,:,:)=Hz(nx-n_cpml_xp+1:nx,:,:)+CPsi_hzx_xp(:,:,:).*Psi_hzx_xp(:,:,:);
% 
% Hx(:,1:n_cpml_yn,:)=Hx(:,1:n_cpml_yn,:)+CPsi_hxy_yn(:,:,:).*Psi_hxy_yn(:,:,:);
% Hz(:,1:n_cpml_yn,:)=Hz(:,1:n_cpml_yn,:)+CPsi_hzy_yn(:,:,:).*Psi_hzy_yn(:,:,:);
% 
% Hx(:,ny-n_cpml_yp+1:ny,:)=Hx(:,ny-n_cpml_yp+1:ny,:)+CPsi_hxy_yp(:,:,:).*Psi_hxy_yp(:,:,:);
% Hz(:,ny-n_cpml_yp+1:ny,:)=Hz(:,ny-n_cpml_yp+1:ny,:)+CPsi_hzy_yp(:,:,:).*Psi_hzy_yp(:,:,:);
% 
% Hx(:,:,1:n_cpml_zn)=Hx(:,:,1:n_cpml_zn)+CPsi_hxz_zn(:,:,:).*Psi_hxz_zn(:,:,:);
% Hy(:,:,1:n_cpml_zn)=Hy(:,:,1:n_cpml_zn)+CPsi_hyz_zn(:,:,:).*Psi_hyz_zn(:,:,:);
% 
% Hx(:,:,nz-n_cpml_zp+1:nz)=Hx(:,:,nz-n_cpml_zp+1:nz)+CPsi_hxz_zp(:,:,:).* Psi_hxz_zp(:,:,:);
% Hy(:,:,nz-n_cpml_zp+1:nz)=Hy(:,:,nz-n_cpml_zp+1:nz)+CPsi_hyz_zp(:,:,:).* Psi_hyz_zp(:,:,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1: n_cpml_xn
Psi_hyx_xn(i,:,:) = cpml_b_mx_xn(i)* Psi_hyx_xn(i,:,:)+cpml_a_mx_xn(i)*...
    ( Ez (i +1,:,:) - Ez(i ,:,:));
Psi_hzx_xn(i,:,:) = cpml_b_mx_xn(i)* Psi_hzx_xn(i,:,:)+cpml_a_mx_xn(i)*...
    ( Ey (i +1,:,:) - Ey(i,:,:));
end
Hy(1:n_cpml_xn,:,:)=Hy(1:n_cpml_xn,:,:)+CPsi_hyx_xn(:,:,:).*Psi_hyx_xn(:,:,:);
Hz(1:n_cpml_xn,:,:)=Hz(1:n_cpml_xn,:,:)+CPsi_hzx_xn(:,:,:).*Psi_hzx_xn(:,:,:);
n_st = nx-n_cpml_xp; 
for i = 1:n_cpml_xp 
Psi_hyx_xp (i,:,:) = cpml_b_mx_xp (i) * Psi_hyx_xp (i ,:,:)+...
    cpml_a_mx_xp (i)*(Ez (i+n_st+1,:,:)-Ez(i+n_st,:,:));
Psi_hzx_xp (i,:,:) = cpml_b_mx_xp (i) * Psi_hzx_xp (i,:,:)+...
    cpml_a_mx_xp (i)*(Ey (i+n_st+1,:,:)-Ey(i+n_st,:,:));
end
Hy(n_st+1:nx,:,:)=Hy(n_st+1:nx,:,:)+CPsi_hyx_xp(:,:,:).*Psi_hyx_xp(:,:,:);
Hz(n_st+1:nx,:,:)=Hz(n_st+1:nx,:,:)+CPsi_hzx_xp(:,:,:).*Psi_hzx_xp(:,:,:);
for j = 1: n_cpml_yn
Psi_hxy_yn(:,j,:)=cpml_b_my_yn(j)*Psi_hxy_yn(:,j,:)+cpml_a_my_yn(j)*...
    (Ez(:,j+1,:) - Ez(:,j,:));
Psi_hzy_yn(:,j,:)=cpml_b_my_yn(j)*Psi_hzy_yn(:,j,:)+cpml_a_my_yn(j)*...
    (Ex(:,j+1,:) - Ex(:,j,:));
end
Hx(:,1:n_cpml_yn,:)=Hx(:,1:n_cpml_yn,:)+CPsi_hxy_yn(:,:,:).*Psi_hxy_yn(:,:,:);
Hz(:,1:n_cpml_yn,:)=Hz(:,1:n_cpml_yn,:)+CPsi_hzy_yn(:,:,:).*Psi_hzy_yn(:,:,:);
n_st = ny-n_cpml_yp; 
for j = 1:n_cpml_yp 
Psi_hxy_yp(:,j,:)=cpml_b_my_yp(j)*Psi_hxy_yp(:,j,:)+cpml_a_my_yp(j)*...
    (Ez(:,j+n_st+1,:)-Ez(:,j+n_st,:));
Psi_hzy_yp(:,j,:)=cpml_b_my_yp(j)*Psi_hzy_yp(:,j,:)+cpml_a_my_yp(j)*...
    (Ex(:,j+n_st+1,:)-Ex(:,j+n_st,:));
end
Hx(:,n_st+1:ny,:)=Hx(:,n_st+1:ny,:)+CPsi_hxy_yp(:,:,:).*Psi_hxy_yp(:,:,:);
Hz(:,n_st+1:ny,:)=Hz(:,n_st+1:ny,:)+CPsi_hzy_yp(:,:,:).*Psi_hzy_yp(:,:,:);
for k = 1: n_cpml_zn
Psi_hxz_zn(:,:,k)=cpml_b_mz_zn(k)*Psi_hxz_zn(:,:,k)+cpml_a_mz_zn(k)*...
    (Ey(:,:,k+1) - Ey(:,:,k));
Psi_hyz_zn(:,:,k)=cpml_b_mz_zn(k)*Psi_hyz_zn(:,:,k)+cpml_a_mz_zn(k)*...
    (Ex(:,:,k+1) - Ex(:,:,k));
end
Hx(:,:,1:n_cpml_zn)=Hx(:,:,1:n_cpml_zn)+CPsi_hxz_zn(:,:,:).*Psi_hxz_zn(:,:,:);
Hy(:,:,1:n_cpml_zn)=Hy(:,:,1:n_cpml_zn)+CPsi_hyz_zn(:,:,:).*Psi_hyz_zn(:,:,:);
n_st = nz-n_cpml_zp; 
for k = 1:n_cpml_zp 
Psi_hyz_zp (:,:,k) = cpml_b_mz_zp(k)*Psi_hyz_zp(:,:,k)+cpml_a_mz_zp(k)*...
    (Ex(:,:,k+n_st+1)-Ex(:,:,k+n_st));
Psi_hxz_zp (:,:,k) = cpml_b_mz_zp(k)*Psi_hxz_zp(:,:,k)+cpml_a_mz_zp(k)*...
    (Ey(:,:,k+n_st+1)-Ey(:,:,k+n_st));
end
Hx(:,:,n_st+1:nz)=Hx(:,:,n_st+1:nz)+CPsi_hxz_zp(:,:,:).* Psi_hxz_zp(:,:,:);
Hy(:,:,n_st+1:nz)=Hy(:,:,n_st+1:nz)+CPsi_hyz_zp(:,:,:).* Psi_hyz_zp(:,:,:);
%% update_magnetic_fields
Hx = Chxh.* Hx+Chxey .*(Ey (1:nxp1,1:ny,2:nzp1)- Ey (1:nxp1,1:ny,1:nz))+...
    Chxez.*(Ez(1:nxp1,2:nyp1,1:nz)- Ez (1:nxp1,1:ny,1:nz));
Hy = Chyh.* Hy+Chyez .*(Ez (2:nxp1,1:nyp1,1:nz)- Ez (1:nx,1:nyp1,1:nz))+...
    Chyex.*(Ex(1:nx,1:nyp1,2:nzp1)- Ex (1:nx,1:nyp1,1:nz));
Hz = Chzh.* Hz+Chzex .*(Ex (1:nx,2:nyp1,1:nzp1)- Ex (1:nx,1:ny,1:nzp1))+...
    Chzey.*(Ey(2:nxp1,1:ny,1:nzp1)- Ey (1:nx,1:ny,1:nzp1));
%% update_electric_fields_for_CPML
% Psi_eyx_xn(:,:,:)= cpml_b_ex_xn.* Psi_eyx_xn(:,:,:)+ cpml_a_ex_xn.*...
%     (Hz(2:n_cpml_xn+1,:,:)-Hz(1:n_cpml_xn,:,:));
% Psi_ezx_xn(:,:,:)= cpml_b_ex_xn.* Psi_ezx_xn(:,:,:)+ cpml_a_ex_xn.*...
%     (Hy(2:n_cpml_xn+1,:,:)-Hy(1:n_cpml_xn,:,:));
% 
% Ey (2:n_cpml_xn+1,:,:) = Ey (2:n_cpml_xn+1,:,:)+CPsi_eyx_xn.*Psi_eyx_xn;
% Ez (2:n_cpml_xn+1,:,:) = Ez (2:n_cpml_xn+1,:,:)+CPsi_ezx_xn.*Psi_ezx_xn; 
% 
% Psi_eyx_xp(:,:,:) = cpml_b_ex_xp.*Psi_eyx_xp(:,:,:)+cpml_a_ex_xp.*...
%     (Hz(1+nx-n_cpml_xp:nx ,:,:) - Hz(nx-n_cpml_xp:nx-1,:,:));
% Psi_ezx_xp(:,:,:) = cpml_b_ex_xp.*Psi_ezx_xp(:,:,:)+cpml_a_ex_xp.*...
%     (Hy(1+nx-n_cpml_xp:nx,:,:) - Hy(nx-n_cpml_xp:nx-1,:,:));
% 
% Ey(nx-n_cpml_xp+1:nx,:,:)=Ey(nx-n_cpml_xp+1:nx,:,:)+CPsi_eyx_xp.*Psi_eyx_xp;
% Ez(nx-n_cpml_xp+1:nx,:,:)=Ez(nx-n_cpml_xp+1:nx,:,:)+CPsi_ezx_xp.*Psi_ezx_xp;
% 
% Psi_exy_yn(:,:,:) = cpml_b_ey_yn.* Psi_exy_yn(:,:,:)+ cpml_a_ey_yn.*...
%     (Hz(:,2:n_cpml_yn+1,:) - Hz(:,1:n_cpml_yn,:));
% Psi_ezy_yn(:,:,:) = cpml_b_ey_yn.* Psi_ezy_yn(:,:,:)+ cpml_a_ey_yn.*...
%     (Hx(:,2:n_cpml_yn+1,:) - Hx(:,1:n_cpml_yn,:));
% 
% Ex (:,2:n_cpml_yn+1,:)=Ex(:,2:n_cpml_yn +1,:)+ CPsi_exy_yn .* Psi_exy_yn;
% Ez (:,2:n_cpml_yn+1,:)=Ez(:,2:n_cpml_yn +1,:)+ CPsi_ezy_yn .* Psi_ezy_yn; 
% 
% Psi_exy_yp(:,:,:) = cpml_b_ey_yp.*Psi_exy_yp(:,:,:)+cpml_a_ey_yp.*...
%     (Hz(:,1+ny-n_cpml_yp:ny,:) - Hz(:,ny-n_cpml_yp:ny-1,:));
% Psi_ezy_yp(:,:,:) = cpml_b_ey_yp.*Psi_ezy_yp(:,:,:)+cpml_a_ey_yp.*...
%     (Hx(:,1+ny-n_cpml_yp:ny,:) - Hx(:,ny-n_cpml_yp:ny-1,:));
% 
% Ex(:,ny-n_cpml_yp+1:ny,:)=Ex(:,ny-n_cpml_yp+1:ny,:)+CPsi_exy_yp.*Psi_exy_yp;
% Ez(:,ny-n_cpml_yp+1:ny,:)=Ez(:,ny-n_cpml_yp+1:ny,:)+CPsi_ezy_yp.*Psi_ezy_yp;
% 
% Psi_exz_zn(:,:,:) = cpml_b_ez_zn.* Psi_exz_zn(:,:,:)+ cpml_a_ez_zn.*...
%     (Hy(:,:,2:n_cpml_zn+1) - Hy(:,:,1:n_cpml_zn));
% Psi_eyz_zn(:,:,:) = cpml_b_ez_zn.* Psi_eyz_zn(:,:,:)+ cpml_a_ez_zn.*...
%     (Hx(:,:,2:n_cpml_zn+1) - Hx(:,:,1:n_cpml_zn));
% 
% Ex (:,:,2:n_cpml_zn+1) = Ex(:,:,2:n_cpml_zn+1)+ CPsi_exz_zn .* Psi_exz_zn;
% Ey (:,:,2:n_cpml_zn +1)= Ey(:,:,2:n_cpml_zn+1)+ CPsi_eyz_zn .* Psi_eyz_zn; 
% 
% Psi_exz_zp(:,:,:) = cpml_b_ez_zp .*Psi_exz_zp(:,:,:)+cpml_a_ez_zp.*...
%     (Hy(:,:,1+nz-n_cpml_zp:nz) - Hy(:,:,nz-n_cpml_zp:nz-1));
% Psi_eyz_zp(:,:,:) = cpml_b_ez_zp.*Psi_eyz_zp(:,:,:)+cpml_a_ez_zp.*...
%     (Hx(:,:,1+nz-n_cpml_zp:nz) - Hx(:,:,nz-n_cpml_zp:nz-1));
% 
% Ex(:,:,nz-n_cpml_zp+1:nz)=Ex(:,:,nz-n_cpml_zp+1:nz)+CPsi_exz_zp.*Psi_exz_zp;
% Ey(:,:,nz-n_cpml_zp+1:nz)=Ey(:,:,nz-n_cpml_zp+1:nz)+CPsi_eyz_zp.*Psi_eyz_zp;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1: n_cpml_xn
Psi_eyx_xn(i,:,:)= cpml_b_ex_xn(i)* Psi_eyx_xn(i,:,:)+ cpml_a_ex_xn(i)*...
    (Hz(i+1,:,:)-Hz(i,:,:));
Psi_ezx_xn(i,:,:)= cpml_b_ex_xn(i)* Psi_ezx_xn(i,:,:)+ cpml_a_ex_xn(i)*...
    (Hy(i+1,:,:)-Hy(i,:,:));
end
Ey (2:n_cpml_xn+1,:,:) = Ey (2:n_cpml_xn+1,:,:)+CPsi_eyx_xn.*Psi_eyx_xn;
Ez (2:n_cpml_xn+1,:,:) = Ez (2:n_cpml_xn+1,:,:)+CPsi_ezx_xn.*Psi_ezx_xn; 
n_st=nx-n_cpml_xp;
for i=1:n_cpml_xp
Psi_eyx_xp(i,:,:) = cpml_b_ex_xp (i)*Psi_eyx_xp(i,:,:)+cpml_a_ex_xp(i)*...
    (Hz(i+n_st ,:,:) - Hz(i+n_st-1,:,:));
Psi_ezx_xp(i,:,:) = cpml_b_ex_xp (i)*Psi_ezx_xp(i,:,:)+cpml_a_ex_xp(i)*...
    (Hy(i+n_st ,:,:) - Hy(i+n_st-1,:,:));
             end
Ey(n_st+1:nx,:,:)=Ey(n_st+1:nx,:,:)+CPsi_eyx_xp.*Psi_eyx_xp;
Ez(n_st+1:nx,:,:)=Ez(n_st+1:nx,:,:)+CPsi_ezx_xp.*Psi_ezx_xp;
for i = 1: n_cpml_yn
Psi_exy_yn(:,i,:) = cpml_b_ey_yn(i) * Psi_exy_yn(:,i,:)+ cpml_a_ey_yn(i)*...
    (Hz(:,i+1,:) - Hz(:,i,:));
Psi_ezy_yn(:,i,:) = cpml_b_ey_yn(i) * Psi_ezy_yn(:,i,:)+ cpml_a_ey_yn(i)*...
    (Hx(:,i+1,:) - Hx(:,i,:));
end
Ex (:,2:n_cpml_yn+1,:)=Ex(:,2:n_cpml_yn +1,:)+ CPsi_exy_yn .* Psi_exy_yn;
Ez (:,2:n_cpml_yn+1,:)=Ez(:,2:n_cpml_yn +1,:)+ CPsi_ezy_yn .* Psi_ezy_yn; 
n_st = ny-n_cpml_yp;
for i=1:n_cpml_yp
Psi_exy_yp(:,i,:) = cpml_b_ey_yp (i)*Psi_exy_yp(:,i,:)+cpml_a_ey_yp(i)*...
    (Hz(:,i+n_st,:) - Hz(:,i+n_st-1,:));
Psi_ezy_yp(:,i,:) = cpml_b_ey_yp (i)*Psi_ezy_yp(:,i,:)+cpml_a_ey_yp(i)*...
    (Hx(:,i+n_st,:) - Hx(:,i+n_st-1,:));
end
Ex(:,n_st+1:ny,:)=Ex(:,n_st+1:ny,:)+CPsi_exy_yp.*Psi_exy_yp;
Ez(:,n_st+1:ny,:)=Ez(:,n_st+1:ny,:)+CPsi_ezy_yp.*Psi_ezy_yp;
for i = 1: n_cpml_zn
Psi_exz_zn(:,:,i) = cpml_b_ez_zn(i) * Psi_exz_zn(:,:,i)+ cpml_a_ez_zn(i)*...
    (Hy(:,:,i+1) - Hy(:,:,i));
Psi_eyz_zn(:,:,i) = cpml_b_ez_zn(i) * Psi_eyz_zn(:,:,i)+ cpml_a_ez_zn(i)*...
    (Hx(:,:,i+1) - Hx(:,:,i));
end
Ex (:,:,2:n_cpml_zn+1) = Ex(:,:,2:n_cpml_zn+1)+ CPsi_exz_zn .* Psi_exz_zn;
Ey (:,:,2:n_cpml_zn +1)= Ey(:,:,2:n_cpml_zn+1)+ CPsi_eyz_zn .* Psi_eyz_zn; 
n_st = nz-n_cpml_zp;
for i =1:n_cpml_zp
Psi_exz_zp(:,:,i) = cpml_b_ez_zp (i)*Psi_exz_zp(:,:,i)+cpml_a_ez_zp(i)*...
    (Hy(:,:,i+n_st) - Hy(:,:,i+n_st-1));
Psi_eyz_zp(:,:,i) = cpml_b_ez_zp (i)*Psi_eyz_zp(:,:,i)+cpml_a_ez_zp(i)*...
    (Hx(:,:,i+n_st) - Hx(:,:,i+n_st-1));
end
Ex(:,:,n_st+1:nz)=Ex(:,:,n_st+1:nz)+CPsi_exz_zp.*Psi_exz_zp;
Ey(:,:,n_st+1:nz)=Ey(:,:,n_st+1:nz)+CPsi_eyz_zp.*Psi_eyz_zp;
%% update_electric_fields
Ex(1:nx,2:ny,2:nz)=Cexe(1:nx,2:ny,2:nz).*Ex(1:nx,2:ny,2:nz)+...
    Cexhz(1:nx,2:ny,2:nz).*(Hz(1:nx,2:ny,2:nz)-Hz(1:nx,1:ny-1,2:nz))+...
    Cexhy(1:nx,2:ny,2:nz).*(Hy(1:nx,2:ny,2:nz)-Hy(1:nx,2:ny,1:nz-1));
Ey(2:nx,1:ny,2:nz)=Ceye(2:nx,1:ny,2:nz).*Ey(2:nx,1:ny,2:nz)+...
    Ceyhx(2:nx,1:ny,2:nz).*(Hx(2:nx,1:ny,2:nz)-Hx(2:nx,1:ny,1:nz-1))+...
    Ceyhz(2:nx,1:ny,2:nz).*(Hz(2:nx,1:ny,2:nz)-Hz(1:nx-1,1:ny,2:nz));
Ez (2:nx,2:ny,1:nz)=Ceze(2:nx,2:ny,1:nz).*Ez(2:nx,2:ny,1:nz)+...
    Cezhy(2:nx,2:ny,1:nz).*(Hy(2:nx,2:ny,1:nz)-Hy(1:nx-1,2:ny,1:nz))+...
    Cezhx(2:nx,2:ny,1:nz).*(Hx(2:nx,2:ny,1:nz)-Hx(2:nx,1:ny-1,1:nz));
%% update source
Ex(1:nx,1:ny,n_cpml_zn+10)=Ex(1:nx,1:ny,n_cpml_zn+10)-sin(omega.*number.*dt);% Differentiated Gaussian pulse
%% capture_and_display_sampled_fields
%     E_0=Ex((nx+1)/2,(ny+1)/2,100);
%     E_1=Ex((nx+1)/2,(ny+1)/2,101);
%     if E_0P~0
%     E_1./E_0P
%     end
% if mod(number,100)==0||number==totalTimeStep
% 
% Hx_sample=Hx(10:end-10,10:end-10,80:end);
% Hy_sample=Hy(10:end-10,10:end-10,80:end);
% Hz_sample=Hz(10:end-10,10:end-10,80:end);
% Ex_sample=Ex(10:end-10,10:end-10,80:end);
% Ey_sample=Ey(10:end-10,10:end-10,80:end);
% Ez_sample=Ez(10:end-10,10:end-10,80:end);

% Hh_1(number/1)=figure


% [X,Z]=meshgrid(xcoor(10:end-10),zcoor(60:end));
% Z=reshape(Z,49,99);
% X=reshape(X,49,99);
% mesh(X,Z,Ex_sample_2D)
% drawnow

% load MyColormap
% colormap(jet)%mymap
% colorbar
% [x,y,z]=meshgrid(xcoor(10:end-10),ycoor(10:ny-10), zcoor(80:end));
% 
% h=slice(x,y,z,Ex_sample,[],[],0)
% % h_p=slice(x,y,z,Ex_sample,[],-7.*10^-3,[])
% h_pp=slice(x,y,z,Ex_sample,-7.*10^-3,[],[])
% set(h, 'edgecolor', 'none')
% % set(h_p, 'edgecolor', 'none')
% set(h_pp, 'edgecolor', 'none')
% caxis([-0.28795 -0.28705])
% drawnow
% caxis([-1 1])
% 
% xlabel('x','fontsize', 12);ylabel('y','fontsize', 12);zlabel('z','fontsize', 12);
% title({'Ex'}); 
% colorbar
% drawnow
% savefig(Hh_1,'three_layer_non_dispersive_skin_fdtd_cpml_gaussian_3D_3D.fig')
% close all
% for i=1:nzp1
%     Ex_sample_2D(i)=Ex(nx/2,ny/2,i) ;
% end
% %  Hh_2=figure
% hold on                                                 
% plot(zcoor(1:end),Ex_sample_2D(1:end))
% xlabel('z','fontsize', 12);ylabel('Ex','fontsize', 12);
% ylim([-1 1]);
% grid on
% drawnow
% savefig(['layers_',num2str(number),'.fig'])
%  end
% disp(number)
end
toc