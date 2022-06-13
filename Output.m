function Output(sNew,p,vx,vy)

global Nf len

x = 1:Nf(1);
y = 1:Nf(2);

figure(01);
DisplayVariable(p,len,1,'Pressure field'), caxis([min(p(:)) max(p(:))]), colorbar;

figure(02);
DisplayVariable(sNew,len,1,'Saturation field'), colorbar;

vex(1:Nf(1),:) = (vx(1:Nf(1),:) + vx(2:Nf(1)+1,:))./2;
vex = (vex)';
vey(:,1:Nf(2)) = (vy(:,1:Nf(2)) + vy(:,2:Nf(2)+1))./2;
vey = (vey)';

figure(05);
quiver(x,y,vex,vey,3),axis([0 Nf(1) 0 Nf(2)]), title('Velocity field'), axis([1 Nf(1) 1 Nf(2)]);

