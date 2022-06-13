function [Disp] = Dispersion(Nf,dx,vx,vy)
%DISPERSION
%-------------------------------------------------------------------------%
%
%              %-----------------------------------------------%
%              %      Rouven Kuenze, University of Lausanne    %
%              %      rouven.kunze@unil.ch                     %
%              %      Ivan Lunati, Univerity of Lausanne       %
%              %      ivan.lunati@unil.ch                      %
%              %-----------------------------------------------%
%
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     PARAMETERS      DEFINITION      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global alphal alphat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Preparing velocities         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qx = vx./dx(2);                                                            % velocity in m/s bcause vx is in m2/s
qy = vy./dx(1);

vyx = sparse(Nf(1)+1,Nf(2));
vxy = sparse(Nf(1),Nf(2)+1);

vyx(2:Nf(1),:) = (qy(1:Nf(1)-1,1:Nf(2)) + qy(1:Nf(1)-1,2:Nf(2)+1) + qy(2:Nf(1),1:Nf(2)) + qy(2:Nf(1),2:Nf(2)+1))./4; % projext y-velocity to vertical boundaries       
vyx(1,:)       = (qy(1,1:Nf(2)) + qy(1,2:Nf(2)+1))./2;
vyx(Nf(1)+1,:) = (qy(Nf(1),1:Nf(2)) + qy(Nf(1),2:Nf(2)+1))./2;

vxy(:,2:Nf(2)) = (qx(1:Nf(1),1:Nf(2)-1) + qx(2:Nf(1)+1,1:Nf(2)-1) + qx(1:Nf(1),2:Nf(2)) + qx(2:Nf(1)+1,2:Nf(2)))./4; % projext x-velocity to horizontal boundaries
vxy(:,1)       = (qx(1:Nf(1),1) + qx(2:Nf(1)+1,1))./2;
vxy(:,Nf(2)+1) = (qx(1:Nf(1),Nf(2)) + qx(2:Nf(1)+1,Nf(2)))./2;

ux = (qx.^2 + vyx.^2).^0.5;                                                % Absolute velocity value x-direction
uy = (vxy.^2 + qy.^2).^0.5;                                                % Absolute velocity value y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Preparing dispersion coefficient   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dxx = alphal.*qx.^2./ux + alphat.*vyx.^2./ux;
Dyy = alphal.*qy.^2./uy + alphat.*vxy.^2./uy;

Dxy             = (alphal - alphat).*qx.*vyx./ux./4;                       % Non diagonal (horizontal boundaries) entries weighted by a fourth
Dxy(:,1)        = Dxy(:,1).*2;                                             % Boundary values only weighted by half
Dxy(:,Nf(2))    = Dxy(:,Nf(2)).*2;
Dxy(1,:)        = 0;
Dxy(Nf(1)+1,:)  = 0;

Dyx             = (alphal - alphat).*vxy.*qy./uy./4;                       % Non diagonal (vertical boundaries) entries weighted by a fourth
Dyx(1,:)        = Dyx(1,:).*2;
Dyx(Nf(1),:)    = Dyx(Nf(1),:).*2;
Dyx(:,1)        = 0;
Dyx(:,Nf(2)+1)  = 0;

Dxx(isnan(Dxx)) = 0;                                                       % Numerical realization for ux = 0 or uy = 0;                                                        
Dyy(isnan(Dyy)) = 0; 
Dxy(isnan(Dxy)) = 0; 
Dyx(isnan(Dyx)) = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Preparing dispersion matrix    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for main diagonal dispersion coefficients %

east  = sparse(Nf(1),Nf(2));
west  = sparse(Nf(1),Nf(2));
north = sparse(Nf(1),Nf(2));
south = sparse(Nf(1),Nf(2));

east (2:Nf(1),1:Nf(2))   = dx(2)*Dxx(2:Nf(1),:)./dx(1); east(1,:)      = 0;     
north(1:Nf(1),2:Nf(2))   = dx(1)*Dyy(:,2:Nf(2))./dx(2); north(:,1)     = 0;
west (1:Nf(1)-1,1:Nf(2)) = dx(2)*Dxx(2:Nf(1),:)./dx(1); west(Nf(1),:)  = 0;
south(1:Nf(1),1:Nf(2)-1) = dx(1)*Dyy(:,2:Nf(2))./dx(2); south(:,Nf(2)) = 0;

Ds      = [south(:) west(:) zeros(prod(Nf),1) east(:) north(:)];
Ds(:,3) = -sum(Ds,2);
Admain   = spdiags(-Ds,[-Nf(1),-1,0,1,Nf(1)],Nf(1)*Nf(2),Nf(1)*Nf(2));

% for Dxy dispersion coefficients %

east = sparse(Nf(1),Nf(2));
west = sparse(Nf(1),Nf(2));

east(1:Nf(1)-1,:) = Dxy(2:Nf(1),:);
west(2:Nf(1),:)   = Dxy(2:Nf(1),:);

Ds    = [-east(:) (-east(:) + west(:)) west(:) east(:) (east(:) - west(:)) -west(:)];
Adxy  = (spdiags(Ds,[-(Nf(1)+1),-Nf(1),-Nf(1)+1,Nf(1)-1,Nf(1),Nf(1)+1],Nf(1)*Nf(2),Nf(1)*Nf(2)))';
 
% for Dyx dispersion coefficients %

lowera   = sparse(Nf(1),Nf(2));                                                                                     
lowerb   = sparse(Nf(1),Nf(2));                                                                                     
uppera   = sparse(Nf(1),Nf(2));
upperb   = sparse(Nf(1),Nf(2));

lowera(1:Nf(1)-1,1:Nf(2)-1) = Dyx(2:Nf(1),  2:Nf(2));
lowerb(2:Nf(1)  ,1:Nf(2)-1) = Dyx(1:Nf(1)-1,2:Nf(2));                                   
uppera(1:Nf(1)-1,2:Nf(2))   = Dyx(2:Nf(1)  ,2:Nf(2));
upperb(2:Nf(1)  ,2:Nf(2))   = Dyx(1:Nf(1)-1,2:Nf(2));

AdyxA   = spdiags([-lowera(:) uppera(:)],[-Nf(1)-1,Nf(1)-1],Nf(1)*Nf(2),Nf(1)*Nf(2));
AdyxB   = spdiags([lowerb(:) -upperb(:)],[-Nf(1)+1,Nf(1)+1],Nf(1)*Nf(2),Nf(1)*Nf(2));

diagL   = sum(AdyxA,2);
diagU   = sum(AdyxB,2);

Adyx   = AdyxA + AdyxB + diag(diagU(1:prod(Nf)-1,1),1) + diag(diagL(2:prod(Nf),1),-1);

% Add values for cells on boundaries %

a = sparse(Nf(1)+1,Nf(2));                                                 % On boundary cells gradient average is build only of 2 gradients. Center values do not vanish and needs to be added                                                
b = sparse(Nf(1),Nf(2)+1);

east  = sparse(Nf(1),Nf(2));
west  = sparse(Nf(1),Nf(2));
north = sparse(Nf(1),Nf(2));
south = sparse(Nf(1),Nf(2));

a(:,1)     = -Dxy(:,1);
a(:,Nf(2)) =  Dxy(:,Nf(2));
b(1,:)     = -Dyx(1,:);
b(Nf(1),:) =  Dyx(Nf(1),:);

east (2:Nf(1),1:Nf(2))   = a(2:Nf(1),:); east(1,:)      = 0;     
north(1:Nf(1),2:Nf(2))   = b(:,2:Nf(2)); north(:,1)     = 0;
west (1:Nf(1)-1,1:Nf(2)) = a(2:Nf(1),:); west(Nf(1),:)  = 0;
south(1:Nf(1),1:Nf(2)-1) = b(:,2:Nf(2)); south(:,Nf(2)) = 0;

Ds       = [south(:) west(:) -east(:) -north(:)];
Adbound  = spdiags(Ds,[-Nf(1),-1,1,Nf(1)],Nf(1)*Nf(2),Nf(1)*Nf(2));
diagonal = sum(Adbound,2);
Adbound  = Adbound + diag(diagonal);

% Create total dispersive matrix %

Disp = Admain + Adxy + Adyx + Adbound;

end
