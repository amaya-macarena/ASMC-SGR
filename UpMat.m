function [Up,rhs] = UpMat(ix,iy,vx,vy,Nf)
%UPMAT
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
% Acknowledgement:  Thanks are due to Hadi Hajibeygi for contributing to 
%                   the very early development of the code. 
%
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     PARAMETERS      DEFINITION      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global FixT Q QT

n    = Nf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTRUCT  DIRECTIONAL ENTRIES  based on  buildPressureSystem.m %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

upw(1:Nf(1),:) = -ix(1:Nf(1),:)   .*vx(1:Nf(1),:);
ups(:,1:Nf(2)) = -iy(:,1:Nf(2))   .*vy(:,1:Nf(2));
upe(1:Nf(1),:) =  (1-ix(2:Nf(1)+1,:)) .*vx(2:Nf(1)+1,:);
upn(:,1:Nf(2)) =  (1-iy(:,2:Nf(2)+1)) .*vy(:,2:Nf(2)+1);

upd(1:Nf(1),1:Nf(2)) = ix(2:Nf(1)+1,1:Nf(2))   .*vx(2:Nf(1)+1,1:Nf(2))...
                      +iy(1:Nf(1),2:Nf(2)+1)   .*vy(1:Nf(1),2:Nf(2)+1)...
                      -(1-ix(1:Nf(1),1:Nf(2))) .*vx(1:Nf(1),1:Nf(2))...
                      -(1-iy(1:Nf(1),1:Nf(2))) .*vy(1:Nf(1),1:Nf(2));

Txeast(2:n(1),:)   = upe(1:n(1)-1,:);    Txeast(1,:)     = 0;
Tynorth(:,2:n(2))  = upn(:,1:n(2)-1);    Tynorth(:,1)    = 0;
Txwest(1:n(1)-1,:) = upw(2:n(1),:);      Txwest(n(1),:)  = 0;
Tysouth(:,1:n(2)-1)= ups(:,2:n(2));      Tysouth(:,n(2)) = 0;

Ds   = [Tysouth(:) Txwest(:) upd(:) Txeast(:) Tynorth(:)];
Up   = spdiags(Ds,[-n(1) -1 0 1 n(1)],n(1)*n(2),n(1)*n(2));

%%%%%%%%%%%%%%%%%%%%%%%%
%%   CONSTRUCT  RHS   %%
%%%%%%%%%%%%%%%%%%%%%%%%

rhs = sparse(zeros(n(1)*n(2),1));

%------------------------------------------------------%
%                   boundary conditions                %
%------------------------------------------------------%

i1   = 1:n(2);                      i2  = n(2) + (1:n(2)); 
i3   = 2*n(2) + (1:n(1));           i4  = 2*n(2) + n(1) + (1:n(1));
ic1  = 1:n(1):prod(n);              ic2 = n(1):n(1):prod(n); 
ic3  = 1:n(1);                      ic4 = (n(2)-1)*n(1)+1:prod(n);

%---------------------%
% Dirichlet Transport %
%---------------------%

iD1  = ic1';              iD2 = ic2'; 
iD3  = ic3';              iD4 = ic4';
t1         = upw(1,:)';   
t2         = upe(n(1),:)';
t3         = ups(:,1);   
t4         = upn(:,n(2));

rhs(iD1,1) = rhs(iD1,1) - t1.*FixT(i1,1);
rhs(iD2,1) = rhs(iD2,1) - t2.*FixT(i2,1);
rhs(iD3,1) = rhs(iD3,1) - t3.*FixT(i3,1);
rhs(iD4,1) = rhs(iD4,1) - t4.*FixT(i4,1);

%%%%%%%%%%%%%%%%%%%%%%%%
%%     source terms   %%
%%%%%%%%%%%%%%%%%%%%%%%%
                                                      
iQ  = sparse((1+sign(Q(:)))/2);                                              % Indicator for in or outflow

out = spdiags((1-iQ).*Q(:),0,prod(n),prod(n));
in  = iQ.*Q(:).*QT(:);

rhs = sparse(rhs + in);
Up  = sparse(Up  - out);
