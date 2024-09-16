function [SysCoor,ElemNode,B1,B2,B3,B4]=get_mesh_filter(L,H,NelL)

%SysCoor = nodes, ElemNode = Topology
%Square grid -L/2<x<L/2, -H/2<y<H/2, 

NelH=floor(NelL/L*H);
X=linspace(-L/2,L/2,NelL+1); Y=linspace(-H/2,H/2,NelH+1); [X100,Y100]=meshgrid(X,Y);
Nel=(size(X100,1)-1)*(size(Y100,2)-1); ElemNode=zeros(Nel,4);  
SysCoor=[X100(:) Y100(:)]; [~,i]=unique(single(SysCoor),'rows'); SysCoor=SysCoor(i,:); el=0;
for i=1:size(X100,1)-1, for j=1:size(Y100,2)-1
   n1=sub2ind(size(X100),i,j); n2=sub2ind(size(X100),i,j+1); n3=sub2ind(size(X100),i+1,j+1); n4=sub2ind(size(X100),i+1,j);   
   el=el+1; ElemNode(el,:)=[n1 n2 n3 n4];
end, end
 
%Get boundaries nodes
%        B2
%     +-------+
%     |       |
% B1  |       | B3
%     +-------+
%        B4
B2=find(SysCoor(:,2)>  H/2-10*eps);   % Nodes at top boundary
B4=find(SysCoor(:,2)< -H/2+10*eps);   % Nodes at bottom boundary
B1=find(SysCoor(:,1)< -L/2+10*eps);   % Nodes at left boundary
B3=find(SysCoor(:,1)>  L/2-10*eps);   % Nodes at right boundary

fprintf(1,'\n\tNo. of elements is %1.0f (%1.0f nodes) \n',size(ElemNode,1),size(SysCoor,1));


