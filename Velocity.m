function [vx,vy] = Velocity(p,Tx,Ty,g)
%VELOCITY Calculates the velocity 
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
% Acknowledgement:  thanks are due to Manav Tyagi and Hadi Hajibeygi for
%                   contributing to the very early development of the code. 
%
%-------------------------------------------------------------------------%


%-------PARAMETERS------%
 
global Fix ibcs 

N = size(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Calculate velocity direction for alpha  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vx = sparse(N(1)+1,N(2));
vy = sparse(N(1),N(2)+1);

%-------------- velocity dirction inside coarse cells --------------------%

vx(2:N(1),:) = -Tx(2:N(1),:).*(p(2:N(1),:)-p(1:N(1)-1,:));                                                       % velocity direction in x direction
vy(:,2:N(2)) = -Ty(:,2:N(2)).*(p(:,2:N(2))-p(:,1:N(2)-1)) - g(:,2:N(2));                                         % velocity direction in y direction

%--------------------- velocity on boundaries ----------------------------%                                 
  
vx(1,:) = -Tx(1,:).*(p(1,:)-Fix(1:N(2))').*ibcs(1:N(2))';                                                        % West Dirichlet
vx(1,:) = vx(1,:)+(Fix(1:N(2)).*~ibcs(1:N(2)))';                                                                 % West Neumann

vx(N(1)+1,:) = -Tx(N(1)+1,:).*(Fix(N(2)+1:2*N(2))'-p(N(1),:)).*ibcs(N(2)+1:2*N(2))';                             % East (adjust Neumann sign!)
vx(N(1)+1,:) = vx(N(1)+1,:)-(Fix(N(2)+1:2*N(2)).*~ibcs(N(2)+1:2*N(2)))';

vy(:,1) = (-Ty(:,1).*(p(:,1)-Fix(2*N(2)+1:2*N(2)+N(1))) - g(:,1)).*ibcs(2*N(2)+1:2*N(2)+N(1));                              % South
vy(:,1) = vy(:,1)+(Fix(2*N(2)+1:2*N(2)+N(1)).*~ibcs(2*N(2)+1:2*N(2)+N(1)));   

vy(:,N(2)+1) = (-Ty(:,N(2)+1).*(Fix(2*N(2)+N(1)+1:2*N(2)+2*N(1))-p(:,N(2))) - g(:,N(2)+1)).*ibcs(2*N(2)+N(1)+1:2*N(2)+2*N(1));   % North (adjust Neumann sign!)
vy(:,N(2)+1) = vy(:,N(2)+1)-(Fix(2*N(2)+N(1)+1:2*N(2)+2*N(1)).*~ibcs(2*N(2)+N(1)+1:2*N(2)+2*N(1)));
