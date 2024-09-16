clear all, close all; drawnow

fprintf(1,'\t******* STARTING DYNAMIC SOLVER OF PLATE WITH HOLE ********* \n\n');

% load quadratic_plate.mat elemnode node B1 B2 B3 B4 


%            B2
%    +---------------+
%    |               |
%    |    +--B6-+    |
% B1 |    |     |    | B3
%    |    B5    B7   |
%    |    +-B8--+    |
%    |               |
%    +---------------+
%           B4
L=1; l=0.2; NelL=40; [node,elemnode,B1,B2,B3,B4,B5,B6,B7,B8]=get_plate;
E=210e9; ny=0.3; rho=7750;

B0=find(node(:,1)==l/2 & node(:,2)==0); %Node at hole x=l/2, y=0
b0x=2*B0-1; b0y=2*B0;

pressure=1e6;

% Constitutive matrix, plane strain
DE=E/(1+ny)/(1-2*ny)*[1-ny  ny   0
                       ny  1-ny  0
                       0   0 (1-2*ny)/2];        

thickness=10;   %Thickness of tube in z-dir
cp=sqrt(E/rho); %Pressure/longitudinal wave speed in the material
t=0;            %Time [s]

c=E/1000; %Damping coefficient classical Ma+Cv+Ku=R [Pa s]
dt=2e-8; %time step sec [s]

name=['Dynamics_plate_c=' num2str(c/E) 'xE'];

fprintf(1,'\n\n\tSolving problem with name %s \n',name); 
fprintf(1,'\tTime step dt=%1.3g µs \n',dt*1e6); 
 
% boundary dofs
b2x=2*B2-1; b2y=2*B2; b4x=2*B4-1; b4y=2*B4; b1x=2*B1-1; b1y=2*B1; b3x=2*B3-1; b3y=2*B3;
b5x=2*B5-1; b5y=2*B5; b6x=2*B6-1; b6y=2*B6; b7x=2*B7-1; b7y=2*B7; b8x=2*B8-1; b8y=2*B8;

%Declare various variables   
ndof=size(node,1)*2; Nel=size(elemnode,1);
un =  zeros(ndof,1);   %displacements
vn =  zeros(ndof,1);   %velocities
an =  zeros(ndof,1);   %accelerations

%Traction boundary conditions because of pressure P
Fext=zeros(ndof,1); bdofs=[b5x; b6y; b7x; b8y]; S=zeros(4,3); S(:,[1 2])=-pressure;
for el = 1:size(elemnode,1)
      n=elemnode(el,:); nn=[2*n(1)-1 2*n(1) 2*n(2)-1 2*n(2) 2*n(3)-1 2*n(3) 2*n(4)-1 2*n(4)];
      [i,j]=intersect(nn,bdofs); if isempty(i)==0; FP=plani4P(node(n,1)',node(n,2)',S,thickness); Fext(i)=Fext(i)+FP(j);end
end
strain=zeros(4,3,Nel);  stress=zeros(4,3,Nel); intpoints=zeros(4,2,Nel);

%Mass matrix M and dampning matrix C
%Use advantage of exact same size of elements
% area Ae of the elements (equal!)
Ae=(L/NelL)^2; Mi=zeros(size(node,1),1); 
for i=1:size(node,1)
    j=find(elemnode(:)==i); 
    Mi(i)=rho*thickness*length(j)*Ae/4;
end
M=sparse(1:ndof); M(1:2:end)=Mi; M(2:2:end)=Mi; M=diag(M); C=M*c/rho;

step=0; 
while 1 %MAIN LOOP
      step=step+1; t=t+dt;
      %1. update displacements
      un=un+dt*vn+0.5*dt^2*an;
      
      %2. update internal forces     
      Fint=zeros(ndof,1);
      for el=1:Nel 
       m=elemnode(el,:); n=[2*m(1)-1 2*m(1) 2*m(2)-1 2*m(2) 2*m(3)-1 2*m(3) 2*m(4)-1 2*m(4)];
       [Fe,Ke,strain(:,:,el),stress(:,:,el),intpoints(:,:,el)]=plane_4_iso(node(m,1)',node(m,2)',DE,un(n)',thickness);
       Fint(n)=Fint(n)+Fe;
      end     
      
      %3. update velocities and accelerations: fully explicit method
      K=M+0.5*dt*C; 
      f=Fext-Fint-C*(vn+0.5*dt*an);
      an1=K\f; 
      vn = vn + dt*0.5*an + dt*0.5*an1;
      an = an1; 
      
      %4. Save some selected variables
      p(step+1)=Fint(b0x)/(L/NelL*thickness);
      u(step+1)=un(b0x);
      tn(step+1)=t;
      dunorm(step+1)=norm(dt*vn + 0.5*dt^2*an)/norm(un);
      
      plot_step=20; % plot at every x step
      if rem(step,plot_step)==0 
         if step==plot_step, fprintf(1,'\nSteps done: '); end
         if rem(step,10*plot_step)==0,fprintf(1,'\n '); end
         fprintf(1,'%1.0f, ',step);

         figure(1);
         seff=sqrt( stress(:,1,:).^2+stress(:,2,:).^2-stress(:,1,:).*stress(:,2,:)+3*stress(:,3,:).^2 ); %effective stress 2D
         subplot(2,2,1),cla
         plot_seff(node,elemnode,intpoints,seff);
         
         subplot(2,2,2),cla, mag=1e5;
         patch('vertices',node+mag*[un(1:2:end) un(2:2:end)],'faces',elemnode,'edgecol','k','facecol',[.8,.9,1]); axis equal;
         axis([-1 1 -1 1]*L/2*1.1); axis equal; axis off; title('deformation (magnified 10^5)','fontsize',14); drawnow
                          
         subplot(2,2,3),cla,
         plot(tn*1e6,p*1e-6,'b-','linewidth',2); 
         k(1)=xlabel('t [µs]'); k(2)=title('pressure on hole edge [MPa]'); set(k,'fontsize',14); ylim([0 1.5]); grid on; drawnow      

         [sx,sy,sxy]=interpolate2nodes(stress,intpoints,node);
         i=find(node(:,1)>0 & abs(node(:,2))<100*eps); [~,j]=sort(node(i,1)); i=i(j);
         subplot(2,2,4),cla,plot(node(i,1),sx(i)*1e-6,'bo-','linewidth',2); 
         k=title(['\sigma_x(x,0) at t=' num2str(1e6*t) ' µs [MPa]']); grid on; k(2)=xlabel('x'); set(k,'fontsize',14); drawnow
      end

      if dunorm(end)<1e-5, break, end
      
end %steps/while

fprintf(1,'\n\tFINISHED. TOTAL %1.0f steps made. \n',step);

eval(['save ' name '_result.mat']);
