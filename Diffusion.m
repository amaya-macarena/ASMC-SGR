function [Di,ri] = Diffusion(Nf,dx,phix,phiy)
%
%MOLECULARDIFFUSION Creates Diffusion Matrix  
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

global Dif ibcD FixT

Dx = dx(2).*phix.*Dif./dx(1);
Dy = dx(1).*phiy.*Dif./dx(2);

Tdeast  = sparse(Nf(1),Nf(2));
Tdwest  = sparse(Nf(1),Nf(2));
Tdsouth = sparse(Nf(1),Nf(2));
Tdnorth = sparse(Nf(1),Nf(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Preparing & Setting of Matrix Di   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tdeast (2:Nf(1),:)   = Dx(2:Nf(1),:); Tdeast(1,:)      = 0;
Tdnorth(:,2:Nf(2))   = Dy(:,2:Nf(2)); Tdnorth(:,1)     = 0;
Tdwest (1:Nf(1)-1,:) = Dx(2:Nf(1),:); Tdwest(Nf(1),:)  = 0;
Tdsouth(:,1:Nf(2)-1) = Dy(:,2:Nf(2)); Tdsouth(:,Nf(2)) = 0;

Ds      = [Tdsouth(:) Tdwest(:) zeros(prod(Nf),1) Tdeast(:) Tdnorth(:)];
Ds(:,3) = -sum(Ds,2);
Di      = spdiags(-Ds,[-Nf(1),-1,0,1,Nf(1)],Nf(1)*Nf(2),Nf(1)*Nf(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Boundary Conditions        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ri   = sparse(prod(Nf),1);

i1   = 1:Nf(2);                      i2  = Nf(2) + (1:Nf(2)); 
i3   = 2*Nf(2) + (1:Nf(1));          i4  = 2*Nf(2) + Nf(1) + (1:Nf(1));

ic1  = 1:Nf(1):prod(Nf);             ic2 = Nf(1):Nf(1):prod(Nf); 
ic3  = 1:Nf(1);                      ic4 = (Nf(2)-1)*Nf(1)+1:prod(Nf);

t1(1:Nf(2),1)  = Dx(1,1:Nf(2))'.*ibcD(i1);   
t2(1:Nf(2),1)  = Dx(Nf(1)+1,1:Nf(2))'.*ibcD(i2);
t3(1:Nf(1),1)  = Dy(1:Nf(1),1) .*ibcD(i3);   
t4(1:Nf(1),1)  = Dy(1:Nf(1),Nf(2)+1) .*ibcD(i4);   

iD  = [ic1';ic2';ic3';ic4'];
tD  = [t1;t2;t3;t4];
Di  = Di + sparse(iD,iD,tD,prod(Nf),prod(Nf));                     

ri(ic1,1) = ri(ic1,1) + t1.*FixT(i1,1);
ri(ic2,1) = ri(ic2,1) + t2.*FixT(i2,1);
ri(ic3,1) = ri(ic3,1) + t3.*FixT(i3,1);
ri(ic4,1) = ri(ic4,1) + t4.*FixT(i4,1);
