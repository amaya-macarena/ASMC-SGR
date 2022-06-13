function [sNew] = Transport(sN,vx,vy,ix,iy,phix,phiy)
%TRANSPORT
%-------------------------------------------------------------------------%
%
%
%              %-----------------------------------------------%
%              %  (c) Rouven Kuenze, University of Lausanne    %
%              %      rouven.kunze@unil.ch                     %
%              %      Ivan Lunati, Univerity of Lausanne       %
%              %      ivan.lunati@unil.ch                      %
%              %-----------------------------------------------%
%
% Acknowledgement:  Thanks are due to Manav Tyagi and Hadi Hajibeygi for
%                   contributing to the very early development of the code. 
%
%-------------------------------------------------------------------------%

global phi dx dt Dif Nf alphal alphat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Build Transport System     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Up,r]  = UpMat(ix,iy,vx,vy,Nf);                                           % Convective Upwind Matrix


if Dif == 0;
   Di = sparse(prod(Nf),prod(Nf));
   ri = sparse(prod(Nf),1);
else
  [Di,ri] = Diffusion(Nf,dx,phix,phiy);                                   % Diffusive Matrix 
end
   
if alphat == 0 && alphal == 0;
   Disp = sparse(prod(Nf),prod(Nf));
else
  [Disp]  = Dispersion(Nf,dx,vx,vy);                                         % Dispersive Matrix
end

Ac      = sparse(phi.*prod(dx)/dt);                                        % Accumulation term

A       = Up + diag(Ac(:)) + Di + Disp;                                    % Stiffness matrix

rhs     = r + ri + sparse(sN(:).*Ac(:));                                   % Right hand side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Solve the system       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sNew=(reshape(A\rhs,Nf));
