function [es,et,eci]=plane_iso4s(ex,ey,D,ed)
% [es,et,eci]=plane_iso4s(ex,ey,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for a 4 node 
%  isoparametric element in plane strain or plane stress. 
%  Coordinates of the Gauss-points are given.
%
% INPUT:   ex = [x1 x2 x3 x4]       element coordinates
%          ey = [y1 y2 y3 y4] 
%                             
%          D                        constitutive matrix (3x3)
%
%          ed = [u1 u2 .. u8]       element displacement vector
%
% OUTPUT: es = [ sigx sigy tauxy ]  element stress matrix
%                                   one row for each Gauss-point
%
%         et = [ epsx epsy gamxy ]  element strain matrix
%                                   one row for each Gauss-point
%
%        eci = [ x  y ]             Gauss-point coordinates in (x,y)
%                                   one row for each Gauss-point
%
%-------------------------------------------------------------

ngp=4; % 2 Gauss-points in each direction (xsi,eta), 2x2=4

%--------- gauss points --------------------------------------
  g1=0.577350269189626; w1=1;
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

   
    es=[]; et=[];

    eci=N*[ex; ey]'; %Coordinates of integration points in (x,y)
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

      ee=B*ed';
      

      
      ss=(D*ee)';
      et=[et; ee'];
      es=[es; ss];

    end 