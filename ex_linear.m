% Example F1 stress/strain/elasticity 2D
clear, close all

% Parameters
t=1e-3;   %thickness
u=0.3e-3; %displacement load

% Constitutive matrix
E=210e9; v=0.3;
D=E/(1-v^2)*[1  v   0  %plane stress
             v  1   0
             0  0 (1-v)/2];

 
%load mesh
%[node,elemnode,B1,B2,B3,B4,BH]=meshrect; 
[node,elemnode,B1,B2,B3,B4,BH]=meshhole;
    
% Define number of freedoms. Here each node have 2 dofs
ndof=2*size(node,1);

% Boundary dofs
b1x=2*B1-1; b1y=2*B1;
b2x=2*B2-1; b2y=2*B2;
b3x=2*B3-1; b3y=2*B3;
b4x=2*B4-1; b4y=2*B4;

% BC
BC=[b1x    0*b1x
    b1y(1) 0
    b3x 0*b3x+u];

% Now we assemble the global stiffness matrix K
Nel=size(elemnode,1);
K=zeros(ndof,ndof);
for el=1:Nel
    ex=node(elemnode(el,:),1)';
    ey=node(elemnode(el,:),2)';
    Ke=plane_iso4ke(ex,ey,t,D); 
    n=elemnode(el,:);
    elemdof=[n(1)*2-1 n(1)*2 n(2)*2-1 n(2)*2 n(3)*2-1 n(3)*2 n(4)*2-1 n(4)*2];   
    K(elemdof,elemdof)=K(elemdof,elemdof)+Ke;
end

% Force load vector
F=zeros(ndof,1);

% Solve!
[U,Q]=solveq(K,F,BC);

L=max(node(:,1)); H=max(node(:,2));

% Compute stresses, strains, integration points
for el=1:Nel
    ex=node(elemnode(el,:),1)';
    ey=node(elemnode(el,:),2)';
    n=elemnode(el,:);
    elemdof=[n(1)*2-1 n(1)*2 n(2)*2-1 n(2)*2 n(3)*2-1 n(3)*2 n(4)*2-1 n(4)*2]; 
    ed=U(elemdof)';
    [stress(:,:,el),strain(:,:,el),intpoints(:,:,el)]=plane_iso4s(ex,ey,D,ed);
end

%Interpolate stresses and strains to nodes (from integration points)
[sx,sy,sxy]=interpolate2nodes(stress,intpoints,node);
[ex,ey,gxy]=interpolate2nodes(strain,intpoints,node);

so=E*u/L; %Theoretical uniaxial stress (Hooke's law)

%post-plotting
%---
figure(1)
subplot(2,1,1)
patch('vertices',node,'faces',elemnode,'facecolor','w','linestyle','-'); 
axis equal; hold on
plot(node(:,1),node(:,2),'ro');
title('Mesh with nodes indicated','fontsize',20)

subplot(2,1,2)
dUx=U(1:2:end); dUy=U(2:2:end); %Nodal displacement
patch('vertices',node,'faces',elemnode,'facecolor','w','linestyle',':'); 
axis equal; hold on
magn=50;
patch('vertices',node+magn*[dUx dUy],'faces',elemnode,'facecolor','b','linestyle','-','facealpha',0.3); 
title(['Deformed (' num2str(magn) 'x magnification) and undeformed meshes'],'fontsize',20)

%---
figure(2)
subplot(3,1,1)
colormap jet
patch('vertices',node,'faces',elemnode,'facevertexcdata',sx,'facecolor','interp','linestyle','none');
figure(2)
axis equal, caxis([-so 4*so]); axis([0 L 0 H]); drawnow; colorbar, title('Stress Sxx [MPa]','fontsize',20); xlabel('x'); ylabel('y');

subplot(3,1,2)
colormap jet
patch('vertices',node,'faces',elemnode,'facevertexcdata',sy,'facecolor','interp','linestyle','none');
figure(2)
axis equal, caxis([-so so]); axis([0 L 0 H]); drawnow; colorbar, title('Stress Syy [MPa]','fontsize',20); xlabel('x'); ylabel('y');

subplot(3,1,3)
colormap jet
patch('vertices',node,'faces',elemnode,'facevertexcdata',sxy,'facecolor','interp','linestyle','none');
figure(2)
axis equal, caxis([-so so]); axis([0 L 0 H]); drawnow; colorbar, title('Stress Sxy [MPa]','fontsize',20); xlabel('x'); ylabel('y'); drawnow

%---
figure(3)

subplot(2,2,1), hold on, grid on, xlabel('x, y=H/2'); ylabel('Sxx');
plot(node(B2,1),sx(B2),'bo');
plot(node(B4,1),sx(B4),'ro');
i=find(node(:,2)==0.01); plot(node(i,1),sx(i),'ko');
ylim([-so 4*so]);
legend('B2','B4','x,y=H/2','Location','best','fontsize',12); legend boxoff

subplot(2,2,2), hold on, grid on, xlabel('x=L/2, y'); ylabel('Sxx');
i=find(node(:,1)==0.03); plot(node(i,2),sx(i),'ks');
ylim([-so 4*so]);
legend('x=L/2,y','fontsize',12,'Location','best'); legend boxoff

subplot(2,2,3), hold on, grid on, xlabel('x=L, y'); ylabel('Sxx');
plot(node(B3,2),sx(B3),'ks'); hold on
ylim([-so 4*so]);
legend('Sxx x=L,y','fontsize',12,'Location','best'); legend boxoff





