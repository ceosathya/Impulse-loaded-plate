function [node,elemnode,B1,B2,B3,B4,BH]=meshrect
% [node,elemnode,B1,B2,B3,B4,BH]=meshrect
%
% node                 : node list (x,y)
% elemnode             : element connectivity matrix [n1 n2 n3 n4]
% B1, B2, B3, B4, BH   : nodes located at respective boundary
%                        BH=[] because not existing
% 
%       +-------------B2-----------+
%       |                          |
%       B1          BH=[]          B3
%       |                          |
%       +-------------B4-----------+
%                   
%   ^ y
%   |
%   +--> x
%
% Per Isaksson 2004-10-07
%

L=0.06; H=0.02;

NelL=20; NelH=NelL*H/L;
X=linspace(0,L,NelL+1); Y=linspace(0,H,NelH+1); [X100,Y100]=meshgrid(X,Y);

Nel=(size(X100,1)-1)*(size(Y100,2)-1); node=[X100(:) Y100(:)]; [~,i]=unique(single(node),'rows'); node=node(i,:); 
elemnode=zeros(Nel,4);  
el=0;
for i=1:size(X100,1)-1
  for j=1:size(Y100,2)-1
   el=el+1;
   n1=sub2ind(size(X100),i,j); n2=sub2ind(size(X100),i,j+1); n3=sub2ind(size(X100),i+1,j+1); n4=sub2ind(size(X100),i+1,j);   
   elemnode(el,:)=[n1 n2 n3 n4];
 end
end

% Get boundaries
%   +----- B2 -----+
%   |              |
%  B1             B3
%   |              |
%   +----- B4 -----+

B1=find( node(:,1)<100*eps ); 
B2=find( node(:,2)>H-100*eps );
B3=find( node(:,1)>L-100*eps );
B4=find( node(:,2)<100*eps );
BH=[];
