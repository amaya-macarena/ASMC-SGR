function DisplayVariable(var,len,scl,stR)

%DISPLAVARIABLE 
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
    
if (nargin<2), len = [1 1]; end
if (nargin<3), scl = 1; end
n = size(var);

x = linspace(0,len(1),2*n(1)+1); x = x(2:2:end);
y = linspace(0,len(2),2*n(2)+1); y = y(2:2:end);
[X,Y] = meshgrid(x,y);

if (size(var,3)==1),
    imagesc(x,y,var')
    set(gca,'YDir','normal');
%     surf(X,Y,var');
else
    quiver(var(:,:,1),var(:,:,2),scl);
end
xlabel('x'); ylabel('y'); 
if (nargin == 4) title(stR); end
view(2);
shading flat
axis image

% plot(var(:,2))