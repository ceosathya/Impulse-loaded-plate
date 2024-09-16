function [Fe,Ke,et,es,eci]=plane_4_iso(ex,ey,C,u,thickness)

ngp=4; g1=1.0/sqrt(3.0); w1=1;
gp(:,1)=[-g1; g1;-g1; g1];  gp(:,2)=[-g1;-g1; g1; g1];
w(:,1)=[ w1; w1; w1; w1];   w(:,2)=[ w1; w1; w1; w1];
wp=w(:,1).*w(:,2); xsi=gp(:,1);  eta=gp(:,2);  r2=ngp*2;

N(:,1)=(1-xsi).*(1-eta)/4;  N(:,2)=(1+xsi).*(1-eta)/4;
N(:,3)=(1+xsi).*(1+eta)/4;  N(:,4)=(1-xsi).*(1+eta)/4;
dNdr(1:2:r2,1)=-(1-eta)/4;     dNdr(1:2:r2,2)= (1-eta)/4;
dNdr(1:2:r2,3)= (1+eta)/4;     dNdr(1:2:r2,4)=-(1+eta)/4;
dNdr(2:2:r2+1,1)=-(1-xsi)/4;   dNdr(2:2:r2+1,2)=-(1+xsi)/4;
dNdr(2:2:r2+1,3)= (1+xsi)/4;   dNdr(2:2:r2+1,4)= (1-xsi)/4;

B=zeros(3,8); Ke=zeros(8,8); Fe=zeros(8,1);

JT=dNdr*[ex;ey]'; et=zeros(4,3); es=zeros(4,3);
eci=N*[ex; ey]'; %Coordinates of integration points in (x,y)

for k=1:4
        indx=[ 2*k-1; 2*k ];
        detJ=det(JT(indx,:)); 
        JTinv=inv(JT(indx,:));
        dNx=JTinv*dNdr(indx,:);
        B(1,1:2:8-1)=dNx(1,:);
        B(2,2:2:8)  =dNx(2,:);
        B(3,1:2:8-1)=dNx(2,:);
        B(3,2:2:8)  =dNx(1,:);        
        de=(B*u')';       
        
        Ke=Ke+B'*C*B*detJ*thickness*wp(k);
        Fe=Fe+B'*C*de'*detJ*thickness*wp(k);      
        
        es(k,:)=(C*de')';  %stress
        et(k,:)=de';       %strain
end


