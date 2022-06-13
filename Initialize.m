function [Tx,Ty,Gr,g] = Initialize(S)
%INITIALIZE
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
global Txk Tyk dx viscosity FixT ibcs gravity density smax
                
N = size(S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 Initialize                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i1 = 1:N(2);             i2  = N(2) + (1:N(2));                            
i3 = 2*N(2) + (1:N(1));  i4  = 2*N(2) + N(1) + (1:N(1));

bcwest(1:N(2))  = FixT(i1);
bceast(1:N(2))  = FixT(i2);
bcsouth(1:N(1)) = FixT(i3);
bcnorth(1:N(1)) = FixT(i4);

Satx           = sparse(N(1)+1,N(2));
Satx(2:N(1),:) = (S(2:N(1),:) + S(1:N(1)-1,:))./2./smax;
Satx(1,:)      = (S(1,:) + bcwest(1:N(2)))./2./smax;
Satx(N(1)+1,:) = (S(N(1),:) + bceast(1:N(2)))./2./smax;

Saty           = sparse(N(1),N(2)+1)./smax;
Saty(:,2:N(2)) = (S(:,2:N(2)) + S(:,1:N(2)-1))./2./smax;
Saty(:,1)      = (S(:,1) + bcsouth(1:N(1))')./2./smax;
Saty(:,N(2)+1) = (S(:,N(2)) + bcnorth(1:N(1))')./2./smax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Calculate upwind Lamda total        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mux   = (1-Satx).*viscosity(1) + Satx.*viscosity(2);
muy   = (1-Saty).*viscosity(1) + Saty.*viscosity(2);

Gy    = ((1-Saty).*density(1) + Saty.*density(2)).*gravity;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   TOTAL MOBILITY and Gravity    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
g            = Tyk./muy.*Gy.*dx(1);                                        % Gravity term with respect to the interfaces
g(:,1)       = g(:,1).*ibcs(2*N(2)+1:2*N(2)+N(1)); 
g(:,N(2)+1)  = g(:,N(2)+1).*ibcs(2*N(2)+N(1)+1:2*N(2)+2*N(1)); 
Gr           = g(:,2:N(2)+1) - g(:,1:N(2));                                % Gravity term corresponding to fine cells

Tx           = Txk./mux.*dx(2)./dx(1);                                     % Transmissivities
Ty           = Tyk./muy.*dx(1)./dx(2);
