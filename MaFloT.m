function [conc]=MaFloT(COND)

%MaFloT Numerical code to solve flow and tracer transport in porous media
%
%  driver(input,output)
%
%  
%  input   ->  must be a string (optional, default 'inputFile')
%  output  ->  must be a string (optional, default 'results/outputFile')
%
%  MaFloT is the main program "driving" the simulation
%  It contains the time loop and calls the following functions:
%  
%  > Initialize
%  > PresMat
%  > Velocity
%  > Transport      >> UpMat
%                   >> Diffusion
%                   >> Dispersion
%
%  > DisplayVariable (is used to display results)
%
%  to see the structure of the InputFile, open InputFile.m
% 
%  see also the following MATLAB functions: reshape, sparse, zeros, ones
%
%%-------------------------------------------------------------------------%
%
%  Conventions used in the code:
%
%  1) fluxes:  outgoing fluxes are negative; ingoing fluxes are positive
%  2) gravity: positive g -> gravity is directed downwards
%  3) scalar fields: p(x,y), s(x,y), Kx(x,y), Ky(x,y), Q(x,y), QT(x,y), ...
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%                             --------------------------
%                             |(1,3)  (2,3) (3,3) (4,3)|
%                          Y  |(1,2)  (2,2) (3,2) (4,2)|
%                             |(1,1)  (2,1) (3,1) (4,1)|
%                             --------------------------
%                                          X
%
%
%  4) boundary conditions (b.c.):
%
%                   ibcs(i) and ibcT(i) define the b.c. type:
%                   1 -> Dirichlet; 0 -> Neumann 
%
%                   Fix(i) and FixT(i) contain the assigned value:
%                   if ibcs(i) = 1 -> Fix(i) is the pressure [Pa]
%                   if ibcs(i) = 0 -> Fix(i) is the flux (per unit 
%                                         of transversal length) [m2/s]
%                   if ibcT(i) = 1 -> FixT(i) is the concentration [kg/m3]
%
%                   the boundary conditions are assigned by a vector that
%                   represents all the cells on the perimeter.
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%
%                              i= 11  i= 12 i= 13 i= 14
%                             --------------------------
%                       i=  3 |(1,3)  (2,3) (3,3) (4,3)| i= 6
%                       i=  2 |(1,2)  (2,2) (3,2) (4,2)| i= 5
%                       i=  1 |(1,1)  (2,1) (3,1) (4,1)| i= 4
%                             --------------------------  
%                              i= 7   i= 8  i=  9 i= 10
%
%  5) vector fields: vx(x,y), vy(x,y), Kx(x,y), Ky(x,y), ...
%                                     [Kx,Ky are diagonal matrix, they can 
%                                              be represented as vector...]
%
%                   vectors are defined at cell interfaces, e.g.,
%                   
%                   vx(i,j) is the velocity between (i,j) and (i+1,j)
%                   vy(i,j) is the velocity between (i,j) and (i,j+1)
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%                   vx:
%                             --------------------------
%                           (1,3)  (2,3) (3,3) (4,3) (5,4)
%                        Y  (1,2)  (2,2) (3,2) (4,2) (5,2)
%                           (1,1)  (2,1) (3,1) (4,1) (5,1)
%                             --------------------------
%                                          X
%                   vy:
%                             |(1,4)  (2,4) (3,4) (4,4)|
%                             |(1,3)  (2,3) (3,3) (4,3)|
%                          Y  |(1,2)  (2,2) (3,2) (4,2)|
%                             |(1,1)  (2,1) (3,1) (4,1)|
%                                          X
%
%  6) linear algebra:
%                   for the solution of the linear system, the scalar
%                   fields are transformed into vectors following the
%                   convention:
%
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%                             --------------------------
%                             |  9     10    11    12  |
%                          Y  |  5      6     7     8  |
%                             |  1      2     3     4  |
%                             --------------------------
%                                          X
%
%                   thus the vector is
%                                       i= 1 |(1,1)|
%                                       i= 2 |(2,1)|
%                                       i= 3 |(3,1)|
%                                       i= 4 |(4,1)|
%                                       i= 5 |(1,2)|
%                                       i= 6 |(2,2)|
%                                       i= 7 |(3,2)|
%                                       i= 8 |(4,2)|
%                                       i= 9 |(1,3)|
%                                       i=10 |(2,3)|
%                                       i=11 |(3,3)|
%                                       i=12 |(4,3)|
%
%
%-------------------------------------------------------------------------%
%
%
%              %-----------------------------------------------%
%              %      Ivan Lunati, Univerity of Lausanne       %
%              %      ivan.lunati@unil.ch                      %
%              %      Rouven Kuenze, University of Lausanne    %
%              %      rouven.kunze@unil.ch                     %
%              %-----------------------------------------------%
%
% Acknowledgement:  thanks to Manav Tyagi, Brad Mallison and Hadi Hajibeygi
%                   for contributing to the early development of the code. 
%
% Maflot is free software; you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation (http://www.gnu.org/licenses/gpl.html);
%-------------------------------------------------------------------------%
%%added by me to force 1 thread computing
LASTN = maxNumCompThreads(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFAULTS AND GLOBAL PARAMETERS  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Nf len Fix Q Txk Tyk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         INITIALIZATION          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InputFile                                                                  % Load Input Data 

% Ela -------------------------------
% K [m2] = COND [m/s] * viscosity/(rho*g) -> 
K = COND * 0.001/(1000*9.81) ;
conv_error='No';
% Ela end ---------------------------

time = dt;                                                    
sNew = s0 .* ones(Nf);

Txk             = sparse(Nf(1)+1,Nf(2));                                   % Interface related permeability (harmonic average)                                                                             
Tyk             = sparse(Nf(1),Nf(2)+1);
Txk(2:Nf(1),:)  = 2./(1./K(1:Nf(1)-1,:) + 1./K(2:Nf(1),:));             
Tyk(:,2:Nf(2))  = 2./(1./K(:,1:Nf(2)-1) + 1./K(:,2:Nf(2)));
Txk(1,:)        = K(1,:);                                    
Txk(Nf(1)+1,:)  = K(Nf(1),:);
Tyk(:,1)        = K(:,1);
Tyk(:,Nf(2)+1)  = K(:,Nf(2));

phix            = sparse(Nf(1)+1,Nf(2));
phiy            = sparse(Nf(1),Nf(2)+1);
phix(2:Nf(1),:) = 2./(1./phi(1:Nf(1)-1,:) + 1./phi(2:Nf(1),:));            % Interface related porosity (harmonic average)
phiy(:,2:Nf(2)) = 2./(1./phi(:,1:Nf(2)-1) + 1./phi(:,2:Nf(2)));
phix(1,:)       = phi(1,:);                                    
phix(Nf(1)+1,:) = phi(Nf(1),:);
phiy(:,1)       = phi(:,1);
phiy(:,Nf(2)+1) = phi(:,Nf(2));

% Ela -------------------------------
conc=zeros(Nf(1),Nf(2),timeSim/dt);
% Ela end ---------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         SIMULATION              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 0;

while (time <= timeSim)                                                    % Time loop
%
  i      = i+1;
  
  % Ela -------------------------------
  %disp(sprintf('Simulation time = %d s',time));
  % Ela end ---------------------------
%
  sOld      = sNew;                                                        % Update saturation of previous iteration (time loop)
  eps       = inf;
  innerIter = 1;
%  
  while (eps >= tolpS && innerIter <= maxpS)                                % Inner loop due to saturation dependence of gravity and viscosity
  %
     sIt = sNew;                                                           % Update saturation of previous iteration (inner loop) 
  %
     [Tx,Ty,G,g]  = Initialize(sNew);
     [A,rhs]      = PresMat(Fix,ibcs,Tx,Ty,G,Q);                           % Construct pressure matrix and rhs vectors
     p            = reshape(A\(rhs),Nf(1),Nf(2));                          % Calculate pressure
  %
     [vx,vy] = Velocity(p,Tx,Ty,g);                                        % Calculate velcities     
  %    
     ix      = (1+sign(vx))/2;                                             % Velocity indicator for upwinding
     iy      = (1+sign(vy))/2;                                             % for aritmetic mean set ix = iy = 1/2
  %
     [sNew]  = Transport(sOld,vx,vy,ix,iy,phix,phiy);                      % Solve transport equation for phase alpha
  %
     eps = norm((abs(sNew(:) - sIt(:))),inf); 
     % Ela -------------------------------
     %disp(sprintf('\t Residual at %d. loop: %d', innerIter,eps));
     % Ela end ---------------------------
     innerIter = innerIter+1;

     % Ela -------------------------------
     %if (innerIter == 40), error('outer finescale loop did not converge'), end
     if (innerIter == 40), disp('outer finescale loop did not converge'),
         conc(:,:,:)=zeros(Nf(1),Nf(2),timeSim/dt)-9999;
         conv_error='Yes';
         time=timeSim;
         break
     end
     % Ela end ---------------------------
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Output             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Ela -------------------------------  
 % Output(sNew,p,vx,vy);
 % Ela end ---------------------------
 % Ela -------------------------------
if strcmp(conv_error,'No');
    conc(:,:,time/dt)=full(sNew);
end
% Ela end ---------------------------
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%     Time Step Control      %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if (time+dt > timeSim) 
     dt = timeSim - time;
     if dt == 0; dt = nan; end
  end
  
  time=time+dt;
  
end