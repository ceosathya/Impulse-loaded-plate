function [Ke]=plane_iso4ke(ex,ey,t,D)
% [Ke]=plane_iso4ke(ex,ey,t,D)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a 4 node isoparametric
%  element in plane strain or plane stress. 
%
% INPUT:  ex = [x1 x2 x3 x4]  element coordinates
%         ey = [y1 y2 y3 y4]
%                             
%         t                   thickness
%
%         D                   constitutive matrix (3x3)
%
%
% OUTPUT: Ke : element stiffness matrix (8 x 8)
%-------------------------------------------------------------

ngp=4; % 2 Gauss-points in each direction (xsi,eta), 2x2=4

%--------- gauss points --------------------------------------
  g1=1/sqrt(3); w1=1;
  gp(:,1)=[-g1; g1;-g1; g1];  gp(:,2)=[-g1;-g1; g1; g1];
  w(:,1)=[ w1; w1; w1; w1];   w(:,2)=[ w1; w1; w1; w1];
  wp=w(:,1).*w(:,2);
  xsi=gp(:,1);  eta=gp(:,2);  r2=ngp*2;

%--------- shape functions -----------------------------------
  N(:,1)=(1-xsi).*(1-eta)/4;  N(:,2)=(1+xsi).*(1-eta)/4;
  N(:,3)=(1+xsi).*(1+eta)/4;  N(:,4)=(1-xsi).*(1+eta)/4;

  dNr(1:2:r2,1)=-(1-eta)/4;     dNr(1:2:r2,2)= (1-eta)/4;
  dNr(1:2:r2,3)= (1+eta)/4;     dNr(1:2:r2,4)=-(1+eta)/4;
  dNr(2:2:r2+1,1)=-(1-xsi)/4;   dNr(2:2:r2+1,2)=-(1+xsi)/4;
  dNr(2:2:r2+1,3)= (1+xsi)/4;   dNr(2:2:r2+1,4)= (1-xsi)/4;

  Ke=zeros(8,8);
  JT=dNr*[ex;ey]';
    
  for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);

      B(1,1:2:8-1)=dNx(1,:);
      B(2,2:2:8)  =dNx(2,:);
      B(3,1:2:8-1)=dNx(2,:);
      B(3,2:2:8)  =dNx(1,:);


      Ke=Ke+B'*D*B*detJ*wp(i)*t;

  end
