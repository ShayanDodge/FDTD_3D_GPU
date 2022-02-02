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

% waveforms.derivative_gaussian.number_of_cells_per_wavelength=20;
% 
% waveforms.cosine_modulated_gaussian(1).bandwith=1e9;
% waveforms.cosine_modulated_gaussian(1).modulation_frequency=waveforms.sinusoidal.frequency;

T=1/waveforms.sinusoidal.frequency;
lambda=(c*T)/1;
omega=2*pi*waveforms.sinusoidal.frequency;
%% FDTD variables
number_of_cells_per_wavelength=(2/7).*150;%2.*50;
dx=lambda/number_of_cells_per_wavelength;
dy=lambda/number_of_cells_per_wavelength;
dz=lambda/number_of_cells_per_wavelength;
totalTime=500*T;
courant_factor=0.86609;%86609
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
boundary.air_buffer_number_of_cells_zp=ceil(1/dz);%20+40+40+40;%+40+200;%+200+200+200+200+200+200+200+200+200;%+200+200+200+200+200+200+200+200
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
%a dielectric_1
material_type(2).eps_r=80;
material_type(2).mu_r=1;
material_type(2).sigma_e=3.3;
material_type(2).sigma_m=0;
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
brick(1).min_x=-0.003;
brick(1).min_y=-0.003;
brick(1).min_z= 0;
brick(1).max_x= 0.003;
brick(1).max_y= 0.003;
brick(1).max_z=2;
brick(1).material_type=1;

% brick(2).min_x=-0.1;
% brick(2).min_y=-0.1;
% brick(2).min_z= 0;
% brick(2).max_x= 0.1;
% brick(2).max_y= 0.1;
% brick(2).max_z= 0.5;
% brick(2).material_type=2;
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
    fdtd_domain.min_y = fdtd_domain.min_y - dy *...
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
nx=(round(fdtd_domain.size_x/dx));
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
%create arrays storing the center coordinates of the cells in
fdtd_domain.cell_center_coordinates_x=single(zeros(nx,ny,nz));
fdtd_domain.cell_center_coordinates_y=single(zeros(nx,ny,nz));
fdtd_domain.cell_center_coordinates_z=single(zeros(nx,ny,nz));
for ind =1:nx
    fdtd_domain.cell_center_coordinates_x(ind,:,:)=...
       ((ind-0.5)*dx+fdtd_domain.min_x);
end
for ind =1:ny
    fdtd_domain.cell_center_coordinates_y(:,ind,:)=...
        (ind-0.5)*dy+fdtd_domain.min_y;
end
for ind =1:nz
    fdtd_domain.cell_center_coordinates_z(:,:,ind)=...
        (ind-0.5)*dz+fdtd_domain.min_z;
end
cx=fdtd_domain.cell_center_coordinates_x;
cy=fdtd_domain.cell_center_coordinates_y;
cz=fdtd_domain.cell_center_coordinates_z;
xcoor=gpuArray(linspace (fdtd_domain. min_x, fdtd_domain.max_x, nxp1));
ycoor=gpuArray(linspace (fdtd_domain. min_y, fdtd_domain.max_y, nyp1));
zcoor=gpuArray(linspace (fdtd_domain. min_z, fdtd_domain.max_z, nzp1));
xcoor=single(xcoor);
ycoor=single(ycoor);
zcoor=single(zcoor);
%% material_3d_space
material_3d_space=single(ones(nx,ny,nz));
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

eps_r_x=single((ones(nx,nyp1,nzp1)));
eps_r_y=single((ones(nxp1,ny,nzp1)));
eps_r_z=single((ones(nxp1,nyp1,nz)));
mu_r_x=single((ones(nxp1,ny,nz)));
mu_r_y=single((ones(nx,nyp1,nz)));
mu_r_z=single((ones(nx,ny,nzp1)));
sigma_e_x=single((zeros(nx,nyp1,nzp1)));
sigma_e_y=single((zeros(nxp1,ny,nzp1)));
sigma_e_z=single((zeros(nxp1,nyp1,nz)));
sigma_m_x=single((zeros(nxp1,ny,nz)));
sigma_m_y=single((zeros(nx,nyp1,nz)));
sigma_m_z=single((zeros(nx,ny,nzp1)));
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
Hx=single((zeros(nxp1,nyp1,nzp1,'gpuArray')));
Hy=single((zeros(nxp1,nyp1,nzp1,'gpuArray')));
Hz=single((zeros(nxp1,nyp1,nzp1,'gpuArray')));
Ex=single((zeros(nxp1,nyp1,nzp1,'gpuArray')));
Ey=single((zeros(nxp1,nyp1,nzp1,'gpuArray')));
Ez=single((zeros(nxp1,nyp1,nzp1,'gpuArray')));
% Exp=single((zeros(nxp1,nyp1,nzp1,'gpuArray')));
% Eyp=single((zeros(nxp1,nyp1,nzp1,'gpuArray')));
% Ezp=single((zeros(nxp1,nyp1,nzp1,'gpuArray')));
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

Psi_eyx_xn = gpuArray(zeros (n_cpml_xn , ny, nzp1));
Psi_ezx_xn = gpuArray(zeros (n_cpml_xn , nyp1, nz));
Psi_hyx_xn = gpuArray(zeros (n_cpml_xn , nyp1, nz));
Psi_hzx_xn = gpuArray(zeros (n_cpml_xn , ny, nzp1 ));

CPsi_eyx_xn = gpuArray(Ceyhz (2:n_cpml_xn +1,:,:)* dx);
CPsi_ezx_xn = gpuArray(Cezhy (2:n_cpml_xn +1,:,:)* dx);
CPsi_hyx_xn = gpuArray(Chyez (1:n_cpml_xn , : , :)*dx);
CPsi_hzx_xn = gpuArray(Chzey (1:n_cpml_xn,:,:) * dx) ;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for i = 1: n_cpml_xn
Ceyhz (i+1,:,:) = Ceyhz (i+1,:,:)/ kappa_ex_xn (i); 
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

Psi_eyx_xp = single(zeros (n_cpml_xp , ny, nzp1));
Psi_ezx_xp = single(zeros (n_cpml_xp , nyp1, nz));
Psi_hyx_xp = single(zeros (n_cpml_xp , nyp1, nz));
Psi_hzx_xp = single(zeros (n_cpml_xp , ny, nzp1 ));

CPsi_eyx_xp = single(gpuArray(Ceyhz (nxp1-n_cpml_xp:nx,:,:))* dx);
CPsi_ezx_xp = single(gpuArray(Cezhy (nxp1-n_cpml_xp:nx,:,:))* dx);
CPsi_hyx_xp = single(gpuArray(Chyez (nxp1-n_cpml_xp:nx,:,:))* dx);
CPsi_hzx_xp = single(gpuArray(Chzey (nxp1-n_cpml_xp:nx,:,:))* dx) ;
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


Psi_exy_yn = single(zeros (nx,n_cpml_yn,nzp1)); 
Psi_ezy_yn = single(zeros (nxp1,n_cpml_yn,nz));
Psi_hxy_yn = single(zeros (nxp1,n_cpml_yn,nz));
Psi_hzy_yn = single(zeros(nx,n_cpml_yn, nzp1 ));
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



Psi_exy_yp = single(zeros (nx,n_cpml_yp,nzp1)); 
Psi_ezy_yp = single(zeros (nxp1,n_cpml_yp,nz));
Psi_hxy_yp = single(zeros (nxp1,n_cpml_yp,nz));
Psi_hzy_yp = single(zeros (nx,n_cpml_yp,nzp1));
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

Psi_eyz_zn = single(zeros (nxp1,ny,n_cpml_zn));
Psi_exz_zn = single(zeros (nx,nyp1,n_cpml_zn));
Psi_hyz_zn = single(zeros (nx,nyp1,n_cpml_zn));
Psi_hxz_zn = single(zeros (nxp1,ny,n_cpml_zn ));
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
Chxey (:,:,j) = (Chxey(:,:,j)/ kappa_mz_zn(j));
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
Psi_eyz_zp = single(zeros (nxp1,ny,n_cpml_zp));
Psi_exz_zp = single(zeros (nx,nyp1,n_cpml_zp));
Psi_hyz_zp = single(zeros (nx,nyp1,n_cpml_zp));
Psi_hxz_zp = single(zeros (nxp1,ny,n_cpml_zp));
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
Chxey(:,:,nz-n_cpml_zp+i) = (Chxey(:,:,nz-n_cpml_zp+i)/kappa_mz_zp(i));
Chyex(:,:,nz-n_cpml_zp+i) = Chyex(:,:,nz-n_cpml_zp+i)/kappa_mz_zp(i);
end



% clear sigma_e_x sigma_e_y sigma_e_z sigma_m_x sigma_m_y sigma_m_z...
%     mu_r_x mu_r_y mu_r_z eps_r_x eps_r_y eps_r_z material_3d_space...
%     fdtd_domain.cell_center_coordinates_x fdtd_domain.cell_center_coordinates_y...
%     fdtd_domain.cell_center_coordinates_z boundary brick material_type fdtd_domain...
%     waveforms alpha_ex_xn alpha_ex_xp alpha_ey_yn alpha_ey_yp alpha_ez_zn alpha_ez_zp...
%     alpha_mx_xn alpha_mx_xp alpha_my_yn alpha_my_yp alpha_mz_zn alpha_mz_zp...
%     kappa_ex_xn kappa_ex_xp kappa_ey_yn kappa_ey_yp kappa_ez_zn kappa_ez_zp...
%     kappa_mx_xn kappa_mx_xp kappa_my_yn kappa_my_yp kappa_mz_zn kappa_mz_zp...
%     rho_e rho_m sigma_ex_xn sigma_ex_xp sigma_ey_yn sigma_ey_yp sigma_ez_zn...
%     sigma_ez_zp sigma_mx_xn sigma_mx_xp sigma_my_yn sigma_my_yp sigma_mz_zn...
%     sigma_mz_zp

%% GPU parameters
n_cpml_zn=abs(boundary.cpml_number_of_cells_zn);
n_cpml_zp=abs(boundary.cpml_number_of_cells_zp);

AA1=single(zeros(nx+1,1,nz+1,'logical','gpuArray'));
AA2=single(zeros(nx+1,ny+1,1,'logical','gpuArray'));
AA4=single(zeros(1,ny+1,nz+1,'logical','gpuArray'));
AA5=single(zeros(nx+1,1,nz,'logical','gpuArray'));
AA6=single(zeros(nx+1,ny+1,1,'logical','gpuArray'));
AA7=single(zeros(nx,ny+1,+1,'logical','gpuArray'));
AA8=single(zeros(1,ny+1,nz+1,'logical','gpuArray'));
AA9=single(zeros(nx,+1,nz+1,'logical','gpuArray'));
AA10=single(zeros(+1,ny+1,nz+1,'logical','gpuArray'));

Cexe=single(cat(1,Cexe,AA4));
Cexhz=single(cat(1,Cexhz,AA4));
Cexhy=single(cat(1,Cexhy,AA4));

Ceye=single(cat(2,Ceye,AA1));
Ceyhx=single(cat(2,Ceyhx,AA1));
Ceyhz=single(cat(2,Ceyhz,AA1));

Ceze=single(cat(3,Ceze,AA2));
Cezhy=single(cat(3,Cezhy,AA2));
Cezhx=single(cat(3,Cezhx,AA2));

Chxh=single(cat(2,Chxh,AA5));
Chxh=single(cat(3,Chxh,AA6));
Chxey=single(cat(2,Chxey,AA5));
Chxey=single(cat(3,Chxey,AA6));
Chxez=single(cat(2,Chxez,AA5));
Chxez=single(cat(3,Chxez,AA6));

Chyh=single(cat(3,Chyh,AA7));
Chyh=single(cat(1,Chyh,AA8));
Chyex=single(cat(3,Chyex,AA7));
Chyex=single(cat(1,Chyex,AA8));
Chyez=single(cat(3,Chyez,AA7));
Chyez=single(cat(1,Chyez,AA8));

Chzh=single(cat(2,Chzh,AA9));
Chzh=single(cat(1,Chzh,AA10));
Chzex=single(cat(2,Chzex,AA9));
Chzex=single(cat(1,Chzex,AA10));
Chzey=single(cat(2,Chzey,AA9));
Chzey=single(cat(1,Chzey,AA10));

Psi_hyx_xn=single(cat(3,Psi_hyx_xn,zeros(20,nyp1,1,'gpuArray')));
Psi_hyx_xn=single(cat(1,Psi_hyx_xn,zeros(nxp1-20,nyp1,nzp1,'gpuArray')));
Psi_hyx_xp=single(cat(3,Psi_hyx_xp,zeros(20,nyp1,1,'gpuArray')));
Psi_hyx_xp=single(cat(1,zeros(nxp1-20,nyp1,nzp1,'gpuArray'),Psi_hyx_xp));
Psi_hzx_xn=single(cat(2,Psi_hzx_xn,zeros(20,1,nzp1,'gpuArray')));
Psi_hzx_xn=single(cat(1,Psi_hzx_xn,zeros(nxp1-20,nyp1,nzp1,'gpuArray')));
Psi_hzx_xp=single(cat(2,Psi_hzx_xp,zeros(20,1,nzp1,'gpuArray')));
Psi_hzx_xp=single(cat(1,zeros(nxp1-20,nyp1,nzp1,'gpuArray'),Psi_hzx_xp));

CPsi_hyx_xn=single(cat(3,CPsi_hyx_xn,zeros(20,nyp1,1,'gpuArray')));
CPsi_hyx_xn=single(cat(1,CPsi_hyx_xn,zeros(nxp1-20,nyp1,nzp1,'gpuArray')));
% CPsi_hyx_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),CPsi_hyx_xn));


CPsi_hyx_xp=single(cat(3,CPsi_hyx_xp,zeros(20,nyp1,1,'gpuArray')));
CPsi_hyx_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),CPsi_hyx_xp));
CPsi_hyx_xp=single(cat(1,CPsi_hyx_xp,zeros(1,nyp1,nzp1,'gpuArray')));

CPsi_hzx_xn=single(cat(2,CPsi_hzx_xn,zeros(20,1,nzp1,'gpuArray')));
CPsi_hzx_xn=single(cat(1,CPsi_hzx_xn,zeros(nxp1-20,nyp1,nzp1,'gpuArray')));
% CPsi_hzx_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),CPsi_hzx_xn));

CPsi_hzx_xp=single(cat(2,CPsi_hzx_xp,zeros(20,1,nzp1,'gpuArray')));
CPsi_hzx_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),CPsi_hzx_xp));
CPsi_hzx_xp=single(cat(1,CPsi_hzx_xp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_b1_mx_xn=single(cpml_b_mx_xn'.*ones(20,ny,nzp1,'gpuArray'));
cpml_b1_mx_xn=single(cat(2,cpml_b1_mx_xn,zeros(20,1,nzp1,'gpuArray')));
cpml_b1_mx_xn=single(cat(1,cpml_b1_mx_xn,zeros(nxp1-20,nyp1,nzp1,'gpuArray')));
% cpml_b1_mx_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),cpml_b1_mx_xn));

cpml_b_mx_xn=single(cpml_b_mx_xn'.*ones(20,nyp1,nz,'gpuArray'));
cpml_b_mx_xn=single(cat(3,cpml_b_mx_xn,zeros(20,nyp1,1,'gpuArray')));
cpml_b_mx_xn=single(cat(1,cpml_b_mx_xn,zeros(nxp1-20,nyp1,nzp1,'gpuArray')));
% cpml_b_mx_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),cpml_b_mx_xn));

cpml_a1_mx_xn=single(cpml_a_mx_xn'.*ones(20,ny,nzp1,'gpuArray'));
cpml_a1_mx_xn=single(cat(2,cpml_a1_mx_xn,zeros(20,1,nzp1,'gpuArray')));
cpml_a1_mx_xn=single(cat(1,cpml_a1_mx_xn,zeros(nxp1-20,nyp1,nzp1,'gpuArray')));
% cpml_a1_mx_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),cpml_a1_mx_xn));

cpml_a_mx_xn=single(cpml_a_mx_xn'.*ones(20,nyp1,nz,'gpuArray'));
cpml_a_mx_xn=single(cat(3,cpml_a_mx_xn,zeros(20,nyp1,1,'gpuArray')));
cpml_a_mx_xn=single(cat(1,cpml_a_mx_xn,zeros(nxp1-20,nyp1,nzp1,'gpuArray')));
% cpml_a_mx_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),cpml_a_mx_xn));

cpml_b1_mx_xp=single(cpml_b_mx_xp'.*ones(20,ny,nzp1,'gpuArray'));
cpml_b1_mx_xp=single(cat(2,cpml_b1_mx_xp,zeros(20,1,nzp1,'gpuArray')));
cpml_b1_mx_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),cpml_b1_mx_xp));
cpml_b1_mx_xp=single(cat(1,cpml_b1_mx_xp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_b_mx_xp=single(cpml_b_mx_xp'.*ones(20,nyp1,nz,'gpuArray'));
cpml_b_mx_xp=single(cat(3,cpml_b_mx_xp,zeros(20,nyp1,1,'gpuArray')));
cpml_b_mx_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),cpml_b_mx_xp));
cpml_b_mx_xp=single(cat(1,cpml_b_mx_xp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a1_mx_xp=single(cpml_a_mx_xp'.*ones(20,ny,nzp1,'gpuArray'));
cpml_a1_mx_xp=single(cat(2,cpml_a1_mx_xp,zeros(20,1,nzp1,'gpuArray')));
cpml_a1_mx_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),cpml_a1_mx_xp));
cpml_a1_mx_xp=single(cat(1,cpml_a1_mx_xp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a_mx_xp=single(cpml_a_mx_xp'.*ones(20,nyp1,nz,'gpuArray'));
cpml_a_mx_xp=single(cat(3,cpml_a_mx_xp,zeros(20,nyp1,1,'gpuArray')));
cpml_a_mx_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),cpml_a_mx_xp));
cpml_a_mx_xp=single(cat(1,cpml_a_mx_xp,zeros(1,nyp1,nzp1,'gpuArray')));

Psi_hxy_yn=single(cat(3,zeros(nxp1,20,1,'gpuArray'),Psi_hxy_yn));
Psi_hxy_yn=single(cat(2,zeros(nxp1,nyp1-20,nzp1,'gpuArray'),Psi_hxy_yn));
Psi_hxy_yp=single(cat(3,Psi_hxy_yp,zeros(nxp1,20,1,'gpuArray')));
Psi_hxy_yp=single(cat(2,Psi_hxy_yp,zeros(nxp1,nyp1-20,nzp1,'gpuArray')));
Psi_hzy_yn=single(cat(1,zeros(1,20,nzp1,'gpuArray'),Psi_hzy_yn));
Psi_hzy_yn=single(cat(2,zeros(nxp1,nyp1-20,nzp1,'gpuArray'),Psi_hzy_yn));
Psi_hzy_yp=single(cat(1,Psi_hzy_yp,zeros(1,20,nzp1,'gpuArray')));
Psi_hzy_yp=single(cat(2,Psi_hzy_yp,zeros(nxp1,nyp1-20,nzp1,'gpuArray')));

CPsi_hxy_yn=single(cat(3,CPsi_hxy_yn,zeros(nxp1,20,1,'gpuArray')));
CPsi_hxy_yn=single(cat(2,CPsi_hxy_yn,zeros(nxp1,nyp1-20,nzp1,'gpuArray')));
% CPsi_hxy_yn=single(cat(2,zeros(nxp1,1,nzp1,'gpuArray'),CPsi_hxy_yn));

CPsi_hxy_yp=single(cat(3,CPsi_hxy_yp,zeros(nxp1,20,1,'gpuArray')));
CPsi_hxy_yp=single(cat(2,zeros(nxp1,ny-20,nzp1,'gpuArray'),CPsi_hxy_yp));
CPsi_hxy_yp=single(cat(2,CPsi_hxy_yp,zeros(nxp1,1,nzp1,'gpuArray')));

CPsi_hzy_yn=single(cat(1,CPsi_hzy_yn,zeros(1,20,nzp1,'gpuArray')));
CPsi_hzy_yn=single(cat(2,CPsi_hzy_yn,zeros(nxp1,nyp1-20,nzp1,'gpuArray')));
% CPsi_hzy_yn=single(cat(2,zeros(nxp1,1,nzp1,'gpuArray'),CPsi_hzy_yn));

CPsi_hzy_yp=single(cat(1,CPsi_hzy_yp,zeros(1,20,nzp1,'gpuArray')));
CPsi_hzy_yp=single(cat(2,zeros(nxp1,ny-20,nzp1,'gpuArray'),CPsi_hzy_yp));
CPsi_hzy_yp=single(cat(2,CPsi_hzy_yp,zeros(nxp1,1,nzp1,'gpuArray')));

cpml_b1_my_yn=single(cpml_b_my_yn.*ones(nx,20,nzp1,'gpuArray'));
cpml_b1_my_yn=single(cat(2,cpml_b1_my_yn,zeros(nx,nyp1-20,nzp1,'gpuArray')));
% cpml_b1_my_yn=single(cat(2,zeros(nx,1,nzp1,'gpuArray'),cpml_b1_my_yn));
cpml_b1_my_yn=single(cat(1,cpml_b1_my_yn,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_b_my_yn=single(cpml_b_my_yn.*ones(nxp1,20,nz,'gpuArray'));
cpml_b_my_yn=single(cat(2,cpml_b_my_yn,zeros(nxp1,nyp1-20,nz,'gpuArray')));
% cpml_b_my_yn=single(cat(2,zeros(nxp1,1,nz,'gpuArray'),cpml_b_my_yn));
cpml_b_my_yn=single(cat(3,cpml_b_my_yn,zeros(nxp1,nyp1,1,'gpuArray')));

cpml_a1_my_yn=single(cpml_a_my_yn.*ones(nx,20,nzp1,'gpuArray'));
cpml_a1_my_yn=single(cat(2,cpml_a1_my_yn,zeros(nx,nyp1-20,nzp1,'gpuArray')));
% cpml_a1_my_yn=single(cat(2,zeros(nx,1,nzp1,'gpuArray'),cpml_a1_my_yn));
cpml_a1_my_yn=single(cat(1,cpml_a1_my_yn,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a_my_yn=single(cpml_a_my_yn.*ones(nxp1,20,nz,'gpuArray'));
cpml_a_my_yn=single(cat(2,cpml_a_my_yn,zeros(nxp1,nyp1-20,nz,'gpuArray')));
% cpml_a_my_yn=single(cat(2,zeros(nxp1,1,nz,'gpuArray'),cpml_a_my_yn));
cpml_a_my_yn=single(cat(3,cpml_a_my_yn,zeros(nxp1,nyp1,1,'gpuArray')));

cpml_b1_my_yp=single(cpml_b_my_yp.*ones(nx,20,nzp1,'gpuArray'));
cpml_b1_my_yp=single(cat(2,zeros(nx,ny-20,nzp1,'gpuArray'),cpml_b1_my_yp));
cpml_b1_my_yp=single(cat(2,cpml_b1_my_yp,zeros(nx,1,nzp1,'gpuArray')));
cpml_b1_my_yp=single(cat(1,cpml_b1_my_yp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_b_my_yp=single(cpml_b_my_yp.*ones(nxp1,20,nz,'gpuArray'));
cpml_b_my_yp=single(cat(2,zeros(nxp1,ny-20,nz,'gpuArray'),cpml_b_my_yp));
cpml_b_my_yp=single(cat(2,cpml_b_my_yp,zeros(nxp1,1,nz,'gpuArray')));
cpml_b_my_yp=single(cat(3,cpml_b_my_yp,zeros(nxp1,nyp1,1,'gpuArray')));

cpml_a1_my_yp=single(cpml_a_my_yp.*ones(nx,20,nzp1,'gpuArray'));
cpml_a1_my_yp=single(cat(2,zeros(nx,ny-20,nzp1,'gpuArray'),cpml_a1_my_yp));
cpml_a1_my_yp=single(cat(2,cpml_a1_my_yp,zeros(nx,1,nzp1,'gpuArray')));
cpml_a1_my_yp=single(cat(1,cpml_a1_my_yp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a_my_yp=single(cpml_a_my_yp.*ones(nxp1,20,nz,'gpuArray'));
cpml_a_my_yp=single(cat(2,zeros(nxp1,ny-20,nz,'gpuArray'),cpml_a_my_yp));
cpml_a_my_yp=single(cat(2,cpml_a_my_yp,zeros(nxp1,1,nz,'gpuArray')));
cpml_a_my_yp=single(cat(3,cpml_a_my_yp,zeros(nxp1,nyp1,1,'gpuArray')));


Psi_hxz_zn=single(cat(2,Psi_hxz_zn,zeros(nxp1,1,n_cpml_zn,'gpuArray')));
Psi_hxz_zn=single(cat(3,Psi_hxz_zn,zeros(nxp1,nyp1,nzp1-n_cpml_zn,'gpuArray')));
Psi_hxz_zp=single(cat(2,Psi_hxz_zp,zeros(nxp1,1,n_cpml_zp,'gpuArray')));
Psi_hxz_zp=single(cat(3,zeros(nxp1,nyp1,nzp1-n_cpml_zp,'gpuArray'),Psi_hxz_zp));
Psi_hyz_zn=single(cat(1,Psi_hyz_zn,zeros(1,nyp1,n_cpml_zn,'gpuArray')));
Psi_hyz_zn=single(cat(3,Psi_hyz_zn,zeros(nxp1,nyp1,nzp1-n_cpml_zn,'gpuArray')));
Psi_hyz_zp=single(cat(1,Psi_hyz_zp,zeros(1,nyp1,n_cpml_zp,'gpuArray')));
Psi_hyz_zp=single(cat(3,zeros(nxp1,nyp1,nzp1-n_cpml_zp,'gpuArray'),Psi_hyz_zp));

CPsi_hxz_zn=single(cat(2,CPsi_hxz_zn,zeros(nxp1,1,n_cpml_zn,'gpuArray')));
CPsi_hxz_zn=single(cat(3,CPsi_hxz_zn,zeros(nxp1,nyp1,nzp1-n_cpml_zn,'gpuArray')));
% CPsi_hxz_zn=single(cat(3,zeros(nxp1,nyp1,1,'gpuArray'),CPsi_hxz_zn));

CPsi_hxz_zp=single(cat(2,CPsi_hxz_zp,zeros(nxp1,1,n_cpml_zp,'gpuArray')));
CPsi_hxz_zp=single(cat(3,zeros(nxp1,nyp1,nz-n_cpml_zp,'gpuArray'),CPsi_hxz_zp));
CPsi_hxz_zp=single(cat(3,CPsi_hxz_zp,zeros(nxp1,nyp1,1,'gpuArray')));

CPsi_hyz_zn=single(cat(1,CPsi_hyz_zn,zeros(1,nyp1,n_cpml_zn,'gpuArray')));
CPsi_hyz_zn=single(cat(3,CPsi_hyz_zn,zeros(nxp1,nyp1,nzp1-n_cpml_zn,'gpuArray')));
% CPsi_hyz_zn=single(cat(3,zeros(nxp1,nyp1,1,'gpuArray'),CPsi_hyz_zn));

CPsi_hyz_zp=single(cat(1,CPsi_hyz_zp,zeros(1,nyp1,n_cpml_zp,'gpuArray')));
CPsi_hyz_zp=single(cat(3,zeros(nxp1,nyp1,nz-n_cpml_zp,'gpuArray'),CPsi_hyz_zp));
CPsi_hyz_zp=single(cat(3,CPsi_hyz_zp,zeros(nxp1,nyp1,1,'gpuArray')));

cpml_b_mz_zn=reshape(cpml_b_mz_zn,1,1,n_cpml_zn);
cpml_b1_mz_zn=single(cpml_b_mz_zn.*ones(nx,nyp1,n_cpml_zn,'gpuArray'));
cpml_b1_mz_zn=single(cat(3,cpml_b1_mz_zn,zeros(nx,nyp1,nzp1-n_cpml_zn,'gpuArray')));
% cpml_b1_mz_zn=single(cat(3,zeros(nx,nyp1,1,'gpuArray'),cpml_b1_mz_zn));
cpml_b1_mz_zn=single(cat(1,cpml_b1_mz_zn,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_b_mz_zn=single(cpml_b_mz_zn.*ones(nxp1,ny,n_cpml_zn,'gpuArray'));
cpml_b_mz_zn=single(cat(3,cpml_b_mz_zn,zeros(nxp1,ny,nzp1-n_cpml_zn,'gpuArray')));
% cpml_b_mz_zn=single(cat(3,zeros(nxp1,ny,1,'gpuArray'),cpml_b_mz_zn));
cpml_b_mz_zn=single(cat(2,cpml_b_mz_zn,zeros(nxp1,1,nzp1,'gpuArray')));

cpml_a_mz_zn=reshape(cpml_a_mz_zn,1,1,n_cpml_zn);
cpml_a1_mz_zn=single(cpml_a_mz_zn.*ones(nx,nyp1,n_cpml_zn,'gpuArray'));
cpml_a1_mz_zn=single(cat(3,cpml_a1_mz_zn,zeros(nx,nyp1,nzp1-n_cpml_zn,'gpuArray')));
% cpml_a1_mz_zn=single(cat(3,zeros(nx,nyp1,1,'gpuArray'),cpml_a1_mz_zn));
cpml_a1_mz_zn=single(cat(1,cpml_a1_mz_zn,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a_mz_zn=single(cpml_a_mz_zn.*ones(nxp1,ny,n_cpml_zn,'gpuArray'));
cpml_a_mz_zn=single(cat(3,cpml_a_mz_zn,zeros(nxp1,ny,nzp1-n_cpml_zn,'gpuArray')));
% cpml_a_mz_zn=single(cat(3,zeros(nxp1,ny,1,'gpuArray'),cpml_a_mz_zn));
cpml_a_mz_zn=single(cat(2,cpml_a_mz_zn,zeros(nxp1,1,nzp1,'gpuArray')));

cpml_b_mz_zp=reshape(cpml_b_mz_zp,1,1,n_cpml_zp);
cpml_b1_mz_zp=single(cpml_b_mz_zp.*ones(nx,nyp1,n_cpml_zp,'gpuArray'));
cpml_b1_mz_zp=single(cat(3,zeros(nx,nyp1,nz-n_cpml_zp,'gpuArray'),cpml_b1_mz_zp));
cpml_b1_mz_zp=single(cat(3,cpml_b1_mz_zp,zeros(nx,nyp1,1,'gpuArray')));
cpml_b1_mz_zp=cat(1,cpml_b1_mz_zp,zeros(1,nyp1,nzp1,'gpuArray'));

cpml_b_mz_zp=single(cpml_b_mz_zp.*ones(nxp1,ny,n_cpml_zp,'gpuArray'));
cpml_b_mz_zp=single(cat(3,zeros(nxp1,ny,nz-n_cpml_zp,'gpuArray'),cpml_b_mz_zp));
cpml_b_mz_zp=single(cat(3,cpml_b_mz_zp,zeros(nxp1,ny,1,'gpuArray')));
cpml_b_mz_zp=single(cat(2,cpml_b_mz_zp,zeros(nxp1,1,nzp1,'gpuArray')));

cpml_a_mz_zp=reshape(cpml_a_mz_zp,1,1,n_cpml_zp);
cpml_a1_mz_zp=single(cpml_a_mz_zp.*ones(nx,nyp1,n_cpml_zp,'gpuArray'));
cpml_a1_mz_zp=single(cat(3,zeros(nx,nyp1,nz-n_cpml_zp,'gpuArray'),cpml_a1_mz_zp));
cpml_a1_mz_zp=single(cat(3,cpml_a1_mz_zp,zeros(nx,nyp1,1,'gpuArray')));
cpml_a1_mz_zp=single(cat(1,cpml_a1_mz_zp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a_mz_zp=single(cpml_a_mz_zp.*ones(nxp1,ny,n_cpml_zp,'gpuArray'));
cpml_a_mz_zp=single(cat(3,zeros(nxp1,ny,nz-n_cpml_zp,'gpuArray'),cpml_a_mz_zp));
cpml_a_mz_zp=single(cat(3,cpml_a_mz_zp,zeros(nxp1,ny,1,'gpuArray')));
cpml_a_mz_zp=single(cat(2,cpml_a_mz_zp,zeros(nxp1,1,nzp1,'gpuArray')));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi_eyx_xn=single(cat(2,Psi_eyx_xn,zeros(20,1,nzp1,'gpuArray')));
Psi_eyx_xn=single(cat(1,Psi_eyx_xn,zeros(nxp1-20,nyp1,nzp1,'gpuArray')));
Psi_eyx_xp=single(cat(2,Psi_eyx_xp,zeros(20,1,nzp1,'gpuArray')));
Psi_eyx_xp=single(cat(1,zeros(nxp1-20,nyp1,nzp1,'gpuArray'),Psi_eyx_xp));
Psi_ezx_xn=single(cat(3,Psi_ezx_xn,zeros(20,nyp1,1,'gpuArray')));
Psi_ezx_xn=single(cat(1,Psi_ezx_xn,zeros(nxp1-20,nyp1,nzp1,'gpuArray')));
Psi_ezx_xp=single(cat(3,Psi_ezx_xp,zeros(20,nyp1,1,'gpuArray')));
Psi_ezx_xp=single(cat(1,zeros(nxp1-20,nyp1,nzp1,'gpuArray'),Psi_ezx_xp));

CPsi_eyx_xn=single(cat(2,CPsi_eyx_xn,zeros(20,1,nzp1,'gpuArray')));
CPsi_eyx_xn=single(cat(1,CPsi_eyx_xn,zeros(nx-20,nyp1,nzp1,'gpuArray')));
CPsi_eyx_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),CPsi_eyx_xn));

CPsi_eyx_xp=single(cat(2,CPsi_eyx_xp,zeros(20,1,nzp1,'gpuArray')));
CPsi_eyx_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),CPsi_eyx_xp));
CPsi_eyx_xp=single(cat(1,CPsi_eyx_xp,zeros(1,nyp1,nzp1,'gpuArray')));

CPsi_ezx_xn=single(cat(3,CPsi_ezx_xn,zeros(20,nyp1,1,'gpuArray')));
CPsi_ezx_xn=single(cat(1,CPsi_ezx_xn,zeros(nx-20,nyp1,nzp1,'gpuArray')));
CPsi_ezx_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),CPsi_ezx_xn));

CPsi_ezx_xp=single(cat(3,CPsi_ezx_xp,zeros(20,nyp1,1,'gpuArray')));
CPsi_ezx_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),CPsi_ezx_xp));
CPsi_ezx_xp=single(cat(1,CPsi_ezx_xp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_b1_ex_xn=cpml_b_ex_xn'.*ones(20,nyp1,nz,'gpuArray');
cpml_b1_ex_xn=single(cat(3,cpml_b1_ex_xn,zeros(20,nyp1,1,'gpuArray')));
cpml_b1_ex_xn=single(cat(1,cpml_b1_ex_xn,zeros(nx-20,nyp1,nzp1,'gpuArray')));
cpml_b1_ex_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),cpml_b1_ex_xn));

cpml_b_ex_xn=single(cpml_b_ex_xn'.*ones(20,ny,nzp1,'gpuArray'));
cpml_b_ex_xn=single(cat(2,cpml_b_ex_xn,zeros(20,1,nzp1,'gpuArray')));
cpml_b_ex_xn=single(cat(1,cpml_b_ex_xn,zeros(nx-20,nyp1,nzp1,'gpuArray')));
cpml_b_ex_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),cpml_b_ex_xn));

cpml_a1_ex_xn=single(cpml_a_ex_xn'.*ones(20,nyp1,nz,'gpuArray'));
cpml_a1_ex_xn=single(cat(3,cpml_a1_ex_xn,zeros(20,nyp1,1,'gpuArray')));
cpml_a1_ex_xn=single(cat(1,cpml_a1_ex_xn,zeros(nx-20,nyp1,nzp1,'gpuArray')));
cpml_a1_ex_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),cpml_a1_ex_xn));

cpml_a_ex_xn=single(cpml_a_ex_xn'.*ones(20,ny,nzp1,'gpuArray'));
cpml_a_ex_xn=single(cat(2,cpml_a_ex_xn,zeros(20,1,nzp1,'gpuArray')));
cpml_a_ex_xn=single(cat(1,cpml_a_ex_xn,zeros(nx-20,nyp1,nzp1,'gpuArray')));
cpml_a_ex_xn=single(cat(1,zeros(1,nyp1,nzp1,'gpuArray'),cpml_a_ex_xn));

cpml_b1_ex_xp=single(cpml_b_ex_xp'.*ones(20,nyp1,nz,'gpuArray'));
cpml_b1_ex_xp=single(cat(3,cpml_b1_ex_xp,zeros(20,nyp1,1,'gpuArray')));
cpml_b1_ex_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),cpml_b1_ex_xp));
cpml_b1_ex_xp=single(cat(1,cpml_b1_ex_xp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_b_ex_xp=single(cpml_b_ex_xp'.*ones(20,ny,nzp1,'gpuArray'));
cpml_b_ex_xp=single(cat(2,cpml_b_ex_xp,zeros(20,1,nzp1,'gpuArray')));
cpml_b_ex_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),cpml_b_ex_xp));
cpml_b_ex_xp=single(cat(1,cpml_b_ex_xp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a1_ex_xp=single(cpml_a_ex_xp'.*ones(20,nyp1,nz,'gpuArray'));
cpml_a1_ex_xp=single(cat(3,cpml_a1_ex_xp,zeros(20,nyp1,1,'gpuArray')));
cpml_a1_ex_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),cpml_a1_ex_xp));
cpml_a1_ex_xp=single(cat(1,cpml_a1_ex_xp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a_ex_xp=single(cpml_a_ex_xp'.*ones(20,ny,nzp1,'gpuArray'));
cpml_a_ex_xp=single(cat(2,cpml_a_ex_xp,zeros(20,1,nzp1,'gpuArray')));
cpml_a_ex_xp=single(cat(1,zeros(nx-20,nyp1,nzp1,'gpuArray'),cpml_a_ex_xp));
cpml_a_ex_xp=single(cat(1,cpml_a_ex_xp,zeros(1,nyp1,nzp1,'gpuArray')));

Psi_exy_yn=single(cat(1,zeros(1,20,nzp1,'gpuArray'),Psi_exy_yn));
Psi_exy_yn=single(cat(2,zeros(nxp1,nyp1-20,nzp1,'gpuArray'),Psi_exy_yn));
Psi_exy_yp=single(cat(1,Psi_exy_yp,zeros(1,20,nzp1,'gpuArray')));
Psi_exy_yp=single(cat(2,Psi_exy_yp,zeros(nxp1,nyp1-20,nzp1,'gpuArray')));
Psi_ezy_yn=single(cat(3,zeros(nxp1,20,1,'gpuArray'),Psi_ezy_yn));
Psi_ezy_yn=single(cat(2,zeros(nxp1,nyp1-20,nzp1,'gpuArray'),Psi_ezy_yn));
Psi_ezy_yp=single(cat(3,Psi_ezy_yp,zeros(nxp1,20,1,'gpuArray')));
Psi_ezy_yp=single(cat(2,Psi_ezy_yp,zeros(nxp1,nyp1-20,nzp1,'gpuArray')));

CPsi_exy_yn=single(cat(1,CPsi_exy_yn,zeros(1,20,nzp1,'gpuArray')));
CPsi_exy_yn=single(cat(2,CPsi_exy_yn,zeros(nxp1,ny-20,nzp1,'gpuArray')));
CPsi_exy_yn=single(cat(2,zeros(nxp1,1,nzp1,'gpuArray'),CPsi_exy_yn));

CPsi_exy_yp=single(cat(1,CPsi_exy_yp,zeros(1,20,nzp1,'gpuArray')));
CPsi_exy_yp=single(cat(2,zeros(nxp1,ny-20,nzp1,'gpuArray'),CPsi_exy_yp));
CPsi_exy_yp=single(cat(2,CPsi_exy_yp,zeros(nxp1,1,nzp1,'gpuArray')));

CPsi_ezy_yn=single(cat(3,CPsi_ezy_yn,zeros(nxp1,20,1,'gpuArray')));
CPsi_ezy_yn=single(cat(2,CPsi_ezy_yn,zeros(nxp1,ny-20,nzp1,'gpuArray')));
CPsi_ezy_yn=single(cat(2,zeros(nxp1,1,nzp1,'gpuArray'),CPsi_ezy_yn));

CPsi_ezy_yp=single(cat(3,CPsi_ezy_yp,zeros(nxp1,20,1,'gpuArray')));
CPsi_ezy_yp=single(cat(2,zeros(nxp1,ny-20,nzp1,'gpuArray'),CPsi_ezy_yp));
CPsi_ezy_yp=single(cat(2,CPsi_ezy_yp,zeros(nxp1,1,nzp1,'gpuArray')));

cpml_b1_ey_yn=single(cpml_b_ey_yn.*ones(nxp1,20,nz,'gpuArray'));
cpml_b1_ey_yn=single(cat(2,cpml_b1_ey_yn,zeros(nxp1,ny-20,nz,'gpuArray')));
cpml_b1_ey_yn=single(cat(2,zeros(nxp1,1,nz,'gpuArray'),cpml_b1_ey_yn));
cpml_b1_ey_yn=single(cat(3,cpml_b1_ey_yn,zeros(nxp1,nyp1,1,'gpuArray')));

cpml_b_ey_yn=single(cpml_b_ey_yn.*ones(nx,20,nzp1,'gpuArray'));
cpml_b_ey_yn=single(cat(2,cpml_b_ey_yn,zeros(nx,ny-20,nzp1,'gpuArray')));
cpml_b_ey_yn=single(cat(2,zeros(nx,1,nzp1,'gpuArray'),cpml_b_ey_yn));
cpml_b_ey_yn=single(cat(1,cpml_b_ey_yn,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a1_ey_yn=single(cpml_a_ey_yn.*ones(nxp1,20,nz,'gpuArray'));
cpml_a1_ey_yn=single(cat(2,cpml_a1_ey_yn,zeros(nxp1,ny-20,nz,'gpuArray')));
cpml_a1_ey_yn=single(cat(2,zeros(nxp1,1,nz,'gpuArray'),cpml_a1_ey_yn));
cpml_a1_ey_yn=single(cat(3,cpml_a1_ey_yn,zeros(nxp1,nyp1,1,'gpuArray')));

cpml_a_ey_yn=single(cpml_a_ey_yn.*ones(nx,20,nzp1,'gpuArray'));
cpml_a_ey_yn=single(cat(2,cpml_a_ey_yn,zeros(nx,ny-20,nzp1,'gpuArray')));
cpml_a_ey_yn=single(cat(2,zeros(nx,1,nzp1,'gpuArray'),cpml_a_ey_yn));
cpml_a_ey_yn=single(cat(1,cpml_a_ey_yn,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_b1_ey_yp=single(cpml_b_ey_yp.*ones(nxp1,20,nz,'gpuArray'));
cpml_b1_ey_yp=single(cat(2,zeros(nxp1,ny-20,nz,'gpuArray'),cpml_b1_ey_yp));
cpml_b1_ey_yp=single(cat(2,cpml_b1_ey_yp,zeros(nxp1,1,nz,'gpuArray')));
cpml_b1_ey_yp=single(cat(3,cpml_b1_ey_yp,zeros(nxp1,nyp1,1,'gpuArray')));

cpml_b_ey_yp=single(cpml_b_ey_yp.*ones(nx,20,nzp1,'gpuArray'));
cpml_b_ey_yp=single(cat(2,zeros(nx,ny-20,nzp1,'gpuArray'),cpml_b_ey_yp));
cpml_b_ey_yp=single(cat(2,cpml_b_ey_yp,zeros(nx,1,nzp1,'gpuArray')));
cpml_b_ey_yp=single(cat(1,cpml_b_ey_yp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a1_ey_yp=single(cpml_a_ey_yp.*ones(nxp1,20,nz,'gpuArray'));
cpml_a1_ey_yp=single(cat(2,zeros(nxp1,ny-20,nz,'gpuArray'),cpml_a1_ey_yp));
cpml_a1_ey_yp=single(cat(2,cpml_a1_ey_yp,zeros(nxp1,1,nz,'gpuArray')));
cpml_a1_ey_yp=single(cat(3,cpml_a1_ey_yp,zeros(nxp1,nyp1,1,'gpuArray')));

cpml_a_ey_yp=single(cpml_a_ey_yp.*ones(nx,20,nzp1,'gpuArray'));
cpml_a_ey_yp=single(cat(2,zeros(nx,ny-20,nzp1,'gpuArray'),cpml_a_ey_yp));
cpml_a_ey_yp=single(cat(2,cpml_a_ey_yp,zeros(nx,1,nzp1,'gpuArray')));
cpml_a_ey_yp=single(cat(1,cpml_a_ey_yp,zeros(1,nyp1,nzp1,'gpuArray')));

Psi_exz_zn=single(cat(1,Psi_exz_zn,zeros(1,nyp1,n_cpml_zn,'gpuArray')));
Psi_exz_zn=single(cat(3,Psi_exz_zn,zeros(nxp1,nyp1,nzp1-n_cpml_zn,'gpuArray')));
Psi_exz_zp=single(cat(1,Psi_exz_zp,zeros(1,nyp1,n_cpml_zp,'gpuArray')));
Psi_exz_zp=single(cat(3,zeros(nxp1,nyp1,nzp1-n_cpml_zp,'gpuArray'),Psi_exz_zp));
Psi_eyz_zn=single(cat(2,Psi_eyz_zn,zeros(nxp1,1,n_cpml_zn,'gpuArray')));
Psi_eyz_zn=single(cat(3,Psi_eyz_zn,zeros(nxp1,nyp1,nzp1-n_cpml_zn,'gpuArray')));
Psi_eyz_zp=single(cat(2,Psi_eyz_zp,zeros(nxp1,1,n_cpml_zp,'gpuArray')));
Psi_eyz_zp=single(cat(3,zeros(nxp1,nyp1,nzp1-n_cpml_zp,'gpuArray'),Psi_eyz_zp));

CPsi_exz_zn=single(cat(1,CPsi_exz_zn,zeros(1,nyp1,n_cpml_zn,'gpuArray')));
CPsi_exz_zn=single(cat(3,CPsi_exz_zn,zeros(nxp1,nyp1,nz-n_cpml_zn,'gpuArray')));
CPsi_exz_zn=single(cat(3,zeros(nxp1,nyp1,1,'gpuArray'),CPsi_exz_zn));

CPsi_exz_zp=single(cat(1,CPsi_exz_zp,zeros(1,nyp1,n_cpml_zp,'gpuArray')));
CPsi_exz_zp=single(cat(3,zeros(nxp1,nyp1,nz-n_cpml_zp,'gpuArray'),CPsi_exz_zp));
CPsi_exz_zp=single(cat(3,CPsi_exz_zp,zeros(nxp1,nyp1,1,'gpuArray')));

CPsi_eyz_zn=single(cat(2,CPsi_eyz_zn,zeros(nxp1,1,n_cpml_zn,'gpuArray')));
CPsi_eyz_zn=single(cat(3,CPsi_eyz_zn,zeros(nxp1,nyp1,nz-n_cpml_zn,'gpuArray')));
CPsi_eyz_zn=single(cat(3,zeros(nxp1,nyp1,1,'gpuArray'),CPsi_eyz_zn));

CPsi_eyz_zp=single(cat(2,CPsi_eyz_zp,zeros(nxp1,1,n_cpml_zp,'gpuArray')));
CPsi_eyz_zp=single(cat(3,zeros(nxp1,nyp1,nz-n_cpml_zp,'gpuArray'),CPsi_eyz_zp));
CPsi_eyz_zp=single(cat(3,CPsi_eyz_zp,zeros(nxp1,nyp1,1,'gpuArray')));

cpml_b_ez_zn=reshape(cpml_b_ez_zn,1,1,n_cpml_zn);

cpml_b1_ez_zn=single(cpml_b_ez_zn.*ones(nxp1,ny,n_cpml_zn,'gpuArray'));
cpml_b1_ez_zn=single(cat(3,cpml_b1_ez_zn,zeros(nxp1,ny,nz-n_cpml_zn,'gpuArray')));
cpml_b1_ez_zn=single(cat(3,zeros(nxp1,ny,1,'gpuArray'),cpml_b1_ez_zn));
cpml_b1_ez_zn=single(cat(2,cpml_b1_ez_zn,zeros(nxp1,1,nzp1,'gpuArray')));


cpml_b_ez_zn=single(cpml_b_ez_zn.*ones(nx,nyp1,n_cpml_zn,'gpuArray'));
cpml_b_ez_zn=single(cat(3,cpml_b_ez_zn,zeros(nx,nyp1,nz-n_cpml_zn,'gpuArray')));
cpml_b_ez_zn=single(cat(3,zeros(nx,nyp1,1,'gpuArray'),cpml_b_ez_zn));
cpml_b_ez_zn=single(cat(1,cpml_b_ez_zn,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a_ez_zn=reshape(cpml_a_ez_zn,1,1,n_cpml_zn);

cpml_a1_ez_zn=single(cpml_a_ez_zn.*ones(nxp1,ny,n_cpml_zn,'gpuArray'));
cpml_a1_ez_zn=single(cat(3,cpml_a1_ez_zn,zeros(nxp1,ny,nz-n_cpml_zn,'gpuArray')));
cpml_a1_ez_zn=single(cat(3,zeros(nxp1,ny,1,'gpuArray'),cpml_a1_ez_zn));
cpml_a1_ez_zn=single(cat(2,cpml_a1_ez_zn,zeros(nxp1,1,nzp1,'gpuArray')));

cpml_a_ez_zn=single(cpml_a_ez_zn.*ones(nx,nyp1,n_cpml_zn,'gpuArray'));
cpml_a_ez_zn=single(cat(3,cpml_a_ez_zn,zeros(nx,nyp1,nz-n_cpml_zn,'gpuArray')));
cpml_a_ez_zn=single(cat(3,zeros(nx,nyp1,1,'gpuArray'),cpml_a_ez_zn));
cpml_a_ez_zn=single(cat(1,cpml_a_ez_zn,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_b_ez_zp=reshape(cpml_b_ez_zp,1,1,n_cpml_zp);

cpml_b1_ez_zp=single(cpml_b_ez_zp.*ones(nxp1,ny,n_cpml_zp,'gpuArray'));
cpml_b1_ez_zp=single(cat(3,zeros(nxp1,ny,nz-n_cpml_zp,'gpuArray'),cpml_b1_ez_zp));
cpml_b1_ez_zp=single(cat(3,cpml_b1_ez_zp,zeros(nxp1,ny,1,'gpuArray')));
cpml_b1_ez_zp=single(cat(2,cpml_b1_ez_zp,zeros(nxp1,1,nzp1,'gpuArray')));

cpml_b_ez_zp=single(cpml_b_ez_zp.*ones(nx,nyp1,n_cpml_zp,'gpuArray'));
cpml_b_ez_zp=single(cat(3,zeros(nx,nyp1,nz-n_cpml_zp,'gpuArray'),cpml_b_ez_zp));
cpml_b_ez_zp=single(cat(3,cpml_b_ez_zp,zeros(nx,nyp1,1,'gpuArray')));
cpml_b_ez_zp=single(cat(1,cpml_b_ez_zp,zeros(1,nyp1,nzp1,'gpuArray')));

cpml_a_ez_zp=reshape(cpml_a_ez_zp,1,1,n_cpml_zp);

cpml_a1_ez_zp=single(cpml_a_ez_zp.*ones(nxp1,ny,n_cpml_zp,'gpuArray'));
cpml_a1_ez_zp=single(cat(3,zeros(nxp1,ny,nz-n_cpml_zp,'gpuArray'),cpml_a1_ez_zp));
cpml_a1_ez_zp=single(cat(3,cpml_a1_ez_zp,zeros(nxp1,ny,1,'gpuArray')));
cpml_a1_ez_zp=single(cat(2,cpml_a1_ez_zp,zeros(nxp1,1,nzp1,'gpuArray')));

cpml_a_ez_zp=single(cpml_a_ez_zp.*ones(nx,nyp1,n_cpml_zp,'gpuArray'));
cpml_a_ez_zp=single(cat(3,zeros(nx,nyp1,nz-n_cpml_zp,'gpuArray'),cpml_a_ez_zp));
cpml_a_ez_zp=single(cat(3,cpml_a_ez_zp,zeros(nx,nyp1,1,'gpuArray')));
cpml_a_ez_zp=single(cat(1,cpml_a_ez_zp,zeros(1,nyp1,nzp1,'gpuArray')));

clear AA5 AA6 AA7 AA8 AA9 AA10

% ex_2=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% ex_3=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% ey_1=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% ey_3=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% ez_1=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% ez_2=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% hx_2=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% hx_3=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% hy_1=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% hy_3=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% hz_1=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% hz_2=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));

SH_1=single(ones(nxp1,nyp1,nzp1,'gpuArray'));
SH_2=single(ones(nxp1,nyp1,nzp1,'gpuArray'));
SH_3=single(ones(nxp1,nyp1,nzp1,'gpuArray'));

SH_1(nxp1,:,:)=0;
% SH_1(:,nyp1,:)=0;
SH_1(:,1,:)=0;
% SH_1(:,:,nzp1)=0;
% SH_1(:,:,1)=0;

SH_2(:,nyp1,:)=0;
% SH_2(nxp1,:,:)=0;
% SH_2(1,:,:)=0;
% SH_2(:,:,nzp1)=0;
% SH_2(:,:,1)=0;

SH_3(:,:,nzp1)=0;
% SH_3(nxp1,:,:)=0;
% SH_3(1,:,:)=0;
% SH_3(:,nyp1,:)=0;
% SH_3(:,1,:)=0;

sh_1=single(ones(nxp1,nyp1,nzp1,'gpuArray'));
sh_2=single(ones(nxp1,nyp1,nzp1,'gpuArray'));
sh_3=single(ones(nxp1,nyp1,nzp1,'gpuArray'));
source=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
source(:,:,n_cpml_zn+ceil(1/dz))=1;

sh_1(:,:,nzp1)=0;
sh_1(:,nyp1,:)=0;

sh_2(nxp1,:,:)=0;
sh_2(:,:,nzp1)=0;
sh_2(:,1,:)=0;


sh_3(:,nyp1,:)=0;
sh_3(nxp1,:,:)=0;
%% initialize_output_parameters_3d
% Hx_sample=single(zeros(nxp1,nyp1,nzp1,'gpuArray'))  ;
% Hy_sample=single(zeros(nxp1,nyp1,nzp1,'gpuArray'))  ;
% Hz_sample=single(zeros(nxp1,nyp1,nzp1,'gpuArray'))  ;
Ex_sample=logical((zeros(nxp1,nyp1,nzp1,'gpuArray')));
% Ey_sample=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% Ez_sample=single(zeros(nxp1,nyp1,nzp1,'gpuArray'));
% Hh_1=figure
% Hh_2=figure
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
% profile on
%% initializing_source_waveforms
time=dt*[0:number_of_time_steps-1].';

% for ind=1:size(waveforms.sinusoidal,2)
%     waveforms.sinusoidal(ind).waveform=...
%     sin(2*pi*waveforms.sinusoidal(ind).frequency*time);
% end
% frequency=waveforms.sinusoidal(1).frequency;
% t_0=0;
% tau=0;

for ind=1:size(waveforms.gaussian,2)
    if waveforms.gaussian(ind).number_of_cells_per_wavelength==0
        nc=number_of_cells_per_wavelength;
    else
    nc=waveforms.gaussian(ind).number_of_cells_per_wavelength;
    end
waveforms.gaussian(ind).maximum_frequency=...
    c/(nc*max([dx,dy,dz]));
tau=(nc*max([dx,dy,dz]))/(2*c);
waveforms.gaussian(ind).tau=tau;
t_0=4.5*waveforms.gaussian(ind).tau;
waveforms.gaussian(ind).t_0=t_0;
waveforms.gaussian(ind).waveform=exp(-((time-t_0)/tau).^2);
end
frequency=waveforms.sinusoidal(1).frequency;
% tau=2e-9;
% t_0=10e-9;

% for ind=1:size(waveforms.derivative_gaussian,2)
%     wfrm=waveforms.derivative_gaussian(ind);
%     if wfrm.number_of_cells_per_wavelength==0
%         nc=number_of_cells_per_wavelength;
%     else
%     nc=waveforms.derivative_gaussian(ind).number_of_cells_per_wavelength;
%     end
% waveforms.derivative_gaussian(ind).maximum_frequency=...
%     c/(nc*max([dx,dy,dz]));
% tau=(nc*max([dx,dy,dz]))/(2*c);
% waveforms.derivative_gaussian(ind).tau=tau;
% t_0=4.5*waveforms.derivative_gaussian(ind).tau;
% waveforms.derivative_gaussian(ind).t_0=t_0;
% waveforms.derivative_gaussian(ind).waveform=...
%     -(sqrt(2*exp(1))/tau)*(time-t_0).*exp(-((time-t_0)/tau).^2);
% end
% frequency=waveforms.sinusoidal(1).frequency;

% for ind=1:size(waveforms.cosine_modulated_gaussian,2)
% frequency=...
% waveforms.cosine_modulated_gaussian(ind).modulation_frequency;
% tau=0.966/waveforms.cosine_modulated_gaussian(ind).bandwith;
% waveforms.cosine_modulated_gaussian(ind).tau=tau;
% t_0=4.5*waveforms.cosine_modulated_gaussian(ind).tau;
% waveforms.cosine_modulated_gaussian(ind).t_0=t_0;
% waveforms.cosine_modulated_gaussian(ind).waveform=...
% cos(2*pi*frequency*(time-t_0)).*exp(-((time-t_0)/tau).^2);
% end

clearvars -except Hx Hy Hz Ex Ey Ez ...
    Psi_hyx_xn cpml_b_mx_xn cpml_b1_mx_xn cpml_a_mx_xn cpml_a1_mx_xn Psi_hzx_xn...
    CPsi_hyx_xn CPsi_hzx_xn Psi_hyx_xp cpml_b_mx_xp cpml_b1_mx_xp cpml_a_mx_xp...
    cpml_a1_mx_xp CPsi_hyx_xp CPsi_hzx_xp Psi_hzx_xp Psi_hxy_yn cpml_b_my_yn ...
    cpml_b1_my_yn cpml_a_my_yn cpml_a1_my_yn Psi_hzy_yn CPsi_hxy_yn CPsi_hzy_yn ...
    Psi_hxy_yp cpml_b_my_yp cpml_b1_my_yp cpml_a_my_yp cpml_a1_my_yp Psi_hzy_yp ...
    CPsi_hxy_yp CPsi_hzy_yp Psi_hxz_zn cpml_b_mz_zn cpml_b1_mz_zn cpml_a_mz_zn ...
    cpml_a1_mz_zn Psi_hyz_zn CPsi_hxz_zn CPsi_hyz_zn Psi_hyz_zp cpml_b_mz_zp ...
    cpml_b1_mz_zp cpml_a_mz_zp cpml_a1_mz_zp Psi_hxz_zp CPsi_hxz_zp CPsi_hyz_zp ...
    Chxh Chxey Chxez Chyh Chyez Chyex Chzh Chzex Chzey dt omega i dy dx...
    sh_1 sh_2 sh_3 AA1 AA2 AA4 Cexe Cexhz Cexhy Ceye Ceyhx Ceyhz Ceze Cezhx Cezhy ...
    Psi_eyx_xn cpml_b_ex_xn cpml_b1_ex_xn cpml_a_ex_xn cpml_a1_ex_xn Psi_ezx_xn CPsi_eyx_xn CPsi_ezx_xn ...
    Psi_eyx_xp cpml_b_ex_xp cpml_b1_ex_xp cpml_a_ex_xp cpml_a1_ex_xp Psi_ezx_xp CPsi_eyx_xp CPsi_ezx_xp ...
    Psi_exy_yn cpml_b_ey_yn cpml_b1_ey_yn cpml_a_ey_yn cpml_a1_ey_yn Psi_ezy_yn CPsi_exy_yn CPsi_ezy_yn ...
    Psi_exy_yp cpml_b_ey_yp cpml_b1_ey_yp cpml_a_ey_yp cpml_a1_ey_yp CPsi_exy_yp CPsi_ezy_yp Psi_ezy_yp ...
    Psi_exz_zn cpml_b_ez_zn cpml_b1_ez_zn cpml_a_ez_zn cpml_a1_ez_zn Psi_eyz_zn CPsi_exz_zn CPsi_eyz_zn ...
    Psi_exz_zp cpml_b_ez_zp cpml_b1_ez_zp cpml_a_ez_zp cpml_a1_ez_zp Psi_eyz_zp CPsi_exz_zp CPsi_eyz_zp ...
    dx dy dz SH_1 SH_2 SH_3 n_cpml_zn nx ny nz zcoor source number_of_time_steps t_0...
    tau frequency

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

% Ex(:,:,n_cpml_zn+61)=Ex(:,:,n_cpml_zn+61)+sin(omega.*number.*dt);

if mod(number,100)==0||number==number_of_time_steps
% close all
% 
% % Hx_sample=Hx;
% % Hy_sample=Hy;
% % Hz_sample=Hz;
% Ex_sample=Ex;
% % Ey_sample=Ey;
% % Ez_sample=Ez;
% 
% Hh_1=figure
% figure(Hh_1);
% % load MyColormap
% colormap(jet)%mymap
% colorbar
% [x,y,z]=meshgrid(xcoor(:),ycoor(:), zcoor(:));
% x=double(x);
% y=double(y);
% z=double(z);
% Ex_sample=double(Ex_sample);
% h=slice(x,y,z,Ex_sample,0,[],[])
% % h_p=slice(x,y,z,Ex_sample,[],-7.*10^-3,[])
% % h_pp=slice(x,y,z,Ex_sample,-7.*10^-3,[],[])
% set(h,'edgecolor', 'none')
% % set(h_p, 'edgecolor', 'none')
% % set(h_pp, 'edgecolor', 'none')
% caxis([-0.28795 -0.28705])
% drawnow
% caxis([-1 1])
% xlabel('x','fontsize', 12);ylabel('y','fontsize', 12);zlabel('z','fontsize', 12);
% title({'Ex'}); 
% colorbar
% drawnow
% print(Hh_1,['layers3d_',num2str(number)],'-djpeg')
% 
for i=1:nz+1
    Ex_sample_2D(i)=Ex(ceil((nx+1)/2),ceil((ny+1)/2),i) ;
end
% % % save(Ex_sample_2D,['Ex_sample_2D',num2str(number)],'-mat')
% % Hh_2=figure
% % figure(Hh_2);
% % hold on                                                 
plot(zcoor(1:end),Ex_sample_2D(1:end))
grid on
% xlabel('z','fontsize', 12);ylabel('Ex','fontsize', 12);
% ylim([0 1]);
% % xlim([-0.012 0]);
% % grid on
drawnow
disp(number)
% print(Hh_2,['layers2d_',num2str(number)],'-dtiffn')
% pause(2)

 end

 end
toc
% delete(par);
% profile report