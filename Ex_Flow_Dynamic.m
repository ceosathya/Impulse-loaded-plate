% Example Heat, Dynamics
clear, close all

% Parameters
k=4; Q=45; To=10; thickness=1;

rho=1;     %density
cp=20;   %heat capacity
t=0;       %time
dt=0.1;   %time step

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
qn=30; dL=L/(length(B2)-1); 
Fb(B2)=-qn*dL; Fb(B2(1))=-qn*dL/2; Fb(B2(end))=-qn*dL/2; %sum(Fb(B2))/L=-30 OK!
Fb(B4)=-qn*dL; Fb(B4(1))=-qn*dL/2; Fb(B4(end))=-qn*dL/2; %sum(Fb(B4))/L=-30 OK!


%Initialize the temperature field
T=zeros(ndof,1); dT=1;  endcrit=1e-6;
% Start time integration!
while norm(dT)/norm(T)>endcrit %Suitable stop criterion

   t=t+dt;

   % Now we assemble the global stiffness matrix K
   K=zeros(ndof,ndof); Fl=zeros(ndof,1);

   % Assemble elements
   for el=1:size(elemnode,1)
     n=elemnode(el,:); 
     ex=node(n,1)'; ey=node(n,2)';
     [ke,ce,fe]=flw2(ex,ey,D,rho,cp,Q,thickness);

     Kt=ke+ce/dt; %Tangential stiffness matrix
     Fe=ce*T(n)/dt + fe ;
     
     K(n,n)=K(n,n)+Kt;
     Fl(n)=Fl(n)+Fe; 
   end

% Now compute global temperature T and flux Aq
T_=T; T=solveq(K,Fl+Fb,BC); dT=T-T_; 

fprintf(1,'  AT STEP %2.0f norm(dT)/norm(T) is %1.2e  \n',round(t/dt),norm(dT)/norm(T));

if rem(round(t/dt),50)==0 % plot T at time t
cla, colormap jet
patch('vertices',node,'faces',elemnode,'facevertexcdata',T,'facecolor','interp','linestyle','none');
axis equal; colorbar, title(['Temperature T [°C] at time ' num2str(t)],'fontsize',20); axis off; caxis([10 20]); drawnow
pause 
end


end

save Ex_Flow_Dynamic.mat %save variables if we want

close all
subplot(2,2,1), cla
colormap jet
patch('vertices',node,'faces',elemnode,'facevertexcdata',T,'facecolor','interp','linestyle','none');
axis equal; colorbar, title('Dynamic solution. T [°C]','fontsize',16); axis off; caxis([10 20]); drawnow

load Ex_Flow_Static.mat T
subplot(2,2,2), cla
colormap jet
patch('vertices',node,'faces',elemnode,'facevertexcdata',T,'facecolor','interp','linestyle','none');
axis equal; colorbar, title('Static solution. T [°C]','fontsize',16); axis off; caxis([10 20]); drawnow

