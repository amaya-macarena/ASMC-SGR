%% GRID PARAMETERS ---------------------------------------------------------------------------------%
global Nf len dx 
len     = [75 101];                      % physical length of the domain in x and y direction [m]
Nf      = [75 101];  	                % number of cells in x and y direction
dx      = len./Nf;                      % cell length [m]
%% PERMEABILITY ------------------------------------------------------------------------------------%
% Ela -------------------------------
% K = ones(Nf(1),Nf(2))*1e-10; 	    % permeability field [m2]
%%% K = COND * viscosity/(rho*g) -> 
% load COND; % log10(m/s)
% K = (10.^COND) * 0.001/(1000*9.81) ;
% Ela end ---------------------------
%% Porosity ----------------------------------------------------------------------------------------%
global phi
phi     = ones(Nf(1),Nf(2))*0.3;	    % porosity field
%% INITIAL CONDITIONS--------------------------------------------------------------------------------%
global s0 smax   
s0   = zeros(Nf(1),Nf(2)) + 0.01;       % Initial saturation (normalized concentration) [-]
smax = 1;                               % maximum concentration [kg/m3] for normalization
% 1 mg/l = 1e-6/1e-3 = 1e-3 kg/m3
%% PHASE PROPERTIES---------------------------------------------------------------------------------%
global viscosity density gravity
viscosity = [0.001 0.001];              % viscosity [kg/m/s] [viscosity(s=0) viscosity(s=smax)]
density   = [1000 1000];                % density [kg/m3] [density(s=0) density(s=smax)]
gravity   = 0; %9.81;                       % gravity acceleration in y [m/s2]
%% SIMULATION PARAMETER FOR TRANSPORT --------------------------------------------------------------%
global dt
timeSim  = 10*86400;  %6*3600/10;%  10*86400                  % total simulation time [s]
dt       = 86400/3;     %2*3600/10; %86400/2                       % time step length [s]
tolpS    = 1.e-5;                     % saturation tolerance on pressure-saturation loop [-]
maxpS    = 50;                          % maximum number of pressure saturation loops to converge
%% BC FLUID ----------------------------------------------------------------------------------------%
global Fix ibcs                                
ibcs = zeros(2*sum(Nf),1);              % type 0:Neumann(N); 1:Dirichlet(D)
Fix  = zeros(2*sum(Nf),1);              % value N [m2/s] (inflow>0); D [Pa]   

%ibcs(1:Nf(2))= 0;% left domain boundary
%ibcs(Nf(2)+1:2*Nf(2)) = 0;% right domain boundary
%ibcs(2*Nf(2)+1:2*Nf(2)+Nf(1)) = 0;% bottom domain boundary
%ibcs(2*Nf(2)+Nf(1)+1:2*Nf(2)+2*Nf(1)) = 0% top domain boundary

ibcs(1:Nf(2)) = 0;
ibcs(Nf(2)+1:2*Nf(2)) = 0;
ibcs(2*Nf(2)+1:2*Nf(2)+Nf(1)) = 1;
ibcs(2*Nf(2)+Nf(1)+1:2*Nf(2)+2*Nf(1)) = 1;

%% BC SOLVENT --------------------------------------------------------------------------------------%
 global FixT
 FixT     = zeros(2*sum(Nf),1);           % normalized concentration of boundary flow [-]

 FixT(2*Nf(2)+[14,30,46,62]) = 1;
 FixT(2*Nf(2)+Nf(1)+[14,30,46,62]) = 1;
%% DIFFUSION AND DISPERSION ------------------------------------------------------------------------%
global Dif ibcD alphal alphat
Dif     = 1e-9;%6.6e-6;                       % [m2/s] molecular diffusion 
ibcD    = zeros(2*sum(Nf),1);           % 1 -> Diffusion on boundary cell
alphal  = 1e-1;                            % longitudinal dispersivity [m]
alphat  = 1e-1;                            % transversal dispersivity [m]
%% SOURCE TERMS ------------------------------------------------------------------------------------%
global Q QT
Q       = zeros(Nf);                    % source term [m2/s]; inflow positive
% Set pumping wells in the center
Q(3:7:73,51)=zeros(1,11)-0.0005;
QT      = zeros(Nf);                    % normalized concentration for source term [-]
