% Example Heat, Static
clear all, close all

% Parameters
k=4; Q=45; To=10; thickness=1;

rho=1;     %density
cp=1000;   %heat capacity

% Constitutive matrix, 2D heat
D=k*[1 0
     0 1];

% Get the square mesh
% nodes(x,y) and topology

% Geometry
%
%
%        B2
%     +-------+
%     |       |
% B1  |       | B3
%     +-------+
%        B4
%    
%    
L=2; NumElemAlongSide=8;
[node,elemnode,B1,B2,B3,B4]=get_mesh_filter(L,L,NumElemAlongSide);
       
% Define number of freedoms. Here each node have 1 dof (scalar valued)
ndof=1*size(node,1);

% Define global load vectors Fb and Fl
Fb=zeros(ndof,1);
Fl=zeros(ndof,1);

% Set BC of T at edge b3
% [dof value]
BC=[ B3 B3*0+To];
    
% Set BC of boundary vector Fb 
qn=30; dL=L/(length(B2)-1); %dL element length
Fb(B2)=-qn*dL; Fb(B2(1))=-qn*dL/2; Fb(B2(end))=-qn*dL/2; %sum(Fb(B2))/L=-30 OK!
Fb(B4)=-qn*dL; Fb(B4(1))=-qn*dL/2; Fb(B4(end))=-qn*dL/2; %sum(Fb(B4))/L=-30 OK!

% Now we assemble the global stiffness matrix K
K=zeros(ndof,ndof); 

% Assemble elements
for el=1:size(elemnode,1)
  n=elemnode(el,:); 
  ex=node(n,1)'; ey=node(n,2)';
  [ke,~,fe]=flw2(ex,ey,D,rho,cp,Q,thickness); %Notice that we get both ke and fe
  K(n,n)=K(n,n)+ke;
  Fl(n)=Fl(n)+fe; 
end

% Now compute global temperature T and flux Aq
[T,R]=solveq(K,Fl+Fb,BC); %T temperature, R boundary flows

% Write out some variables
if 1==2
T
R
K
Fb
Fl
end

save Ex_Flow_Static.mat T %save variables if we want

% post plot
colormap jet
patch('vertices',node,'faces',elemnode,'facevertexcdata',T,'facecolor','interp','linestyle','none');
axis equal; colorbar, title('Static solution. Temperature T [°C]','fontsize',20); axis off; caxis([10 20]);


