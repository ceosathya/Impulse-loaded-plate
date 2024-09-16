function [node,elemnode,B1,B2,B3,B4,BH]=meshhole
% [node,elemnode,B1,B2,B3,B4,BH]=meshhole
%
% node                 : node list (x,y)
% elemnode             : element connectivity matrix [n1 n2 n3 n4]
% B1, B2, B3, B4, BH   : nodes located at respective boundary
%                   
%          +-------B2-------+
%          |       --       |
%          |    /      BH   |
%         B1   |        |   B3
%          |    \      /    |
%          |       --       |
%   ^ x2   +-------B4-------+
%   |
%   +--> x1
%

%First make a symmetric 1/4 part
a=0.005; %radius
Nno=65;
Nel=48;

% parametric definition of node coordinates
ba = 1.25*a;
b  = 1.5*a;
bb = 1.75*a;
c  = 2*a;
node = [ 
      a                      0
      cos(pi/16)*a           sin(pi/16)*a
      cos(pi/8)*a            sin(pi/8)*a
      cos(3*pi/16)*a         sin(3*pi/16)*a
      cos(pi/4)*a            sin(pi/4)*a
      cos(5*pi/16)*a         sin(5*pi/16)*a
      cos(3*pi/8)*a          sin(3*pi/8)*a
      cos(7*pi/16)*a         sin(7*pi/16)*a
      0                      a
      ba                     0
      cos(pi/16)*ba          sin(pi/16)*ba
      cos(pi/8)*ba*1.025     sin(pi/8)*ba*1.025
      cos(3*pi/16)*ba*1.05   sin(3*pi/16)*ba*1.05
      cos(pi/4)*ba*1.075     sin(pi/4)*ba*1.075
      cos(5*pi/16)*ba*1.05   sin(5*pi/16)*ba*1.05
      cos(3*pi/8)*ba*1.025   sin(3*pi/8)*ba*1.025
      cos(7*pi/16)*ba        sin(7*pi/16)*ba
      0                      ba
      b                      0
      cos(pi/16)*b           sin(pi/16)*b
      cos(pi/8)*b*1.05       sin(pi/8)*b*1.05
      cos(3*pi/16)*b*1.1     sin(3*pi/16)*b*1.1
      cos(pi/4)*b*1.175      sin(pi/4)*b*1.175
      cos(5*pi/16)*b*1.1     sin(5*pi/16)*b*1.1
      cos(3*pi/8)*b*1.05     sin(3*pi/8)*b*1.05
      cos(7*pi/16)*b         sin(7*pi/16)*b
      0                      b
      bb                     0
      cos(pi/16)*bb          sin(pi/16)*bb
      cos(pi/8)*bb*1.05      sin(pi/8)*bb*1.05
      cos(3*pi/16)*bb*1.15   sin(3*pi/16)*bb*1.15
      cos(pi/4)*bb*1.275     sin(pi/4)*bb*1.275
      cos(5*pi/16)*bb*1.15   sin(5*pi/16)*bb*1.15
      cos(3*pi/8)*bb*1.05    sin(3*pi/8)*bb*1.05
      cos(7*pi/16)*bb        sin(7*pi/16)*bb
      0                      bb
      c                      0
      c                      a/2*0.8
      c                      a*0.85
      c                      3*a/2*0.9
      c                      c
      3*a/2*0.9              c
      a*0.85                 c
      a/2*0.8                c
      0                      c 
      3*a                    0
      3*a                    a/2         
      3*a                    a  
      3*a                    3*a/2       
      3*a                    2*a         
      4*a                    0
      4*a                    a/2         
      4*a                    a  
      4*a                    3*a/2       
      4*a                    2*a         
      5*a                    0
      5*a                    a/2         
      5*a                    a  
      5*a                    3*a/2       
      5*a                    2*a         
      6*a                    0
      6*a                    a/2         
      6*a                    a  
      6*a                    3*a/2       
      6*a                    2*a          ];

% topology 
elemnode = [1 10 11  2;
     2 11 12  3;
     3 12 13  4;
     4 13 14  5;
       5 14 15  6;
       6 15 16  7;
       7 16 17  8;
       8 17 18  9;
      10 19 20 11;
     11 20 21 12;
     12 21 22 13;
     13 22 23 14;    
     14 23 24 15;    
     15 24 25 16;    
     16 25 26 17;    
     17 26 27 18;    
     19 28 29 20;   
     20 29 30 21;   
     21 30 31 22;    
     22 31 32 23;    
     23 32 33 24;    
     24 33 34 25;    
     25 34 35 26;    
     26 35 36 27;    
     28 37 38 29;
     29 38 39 30;
     30 39 40 31;
     31 40 41 32; 
     32 41 42 33;
     33 42 43 34;
     34 43 44 35;
     35 44 45 36;
     37 46 47 38;
     38 47 48 39;
     39 48 49 40;
     40 49 50 41;
     46 51 52 47;
     47 52 53 48;
     48 53 54 49;
     49 54 55 50;
     51 56 57 52;
     52 57 58 53;
     53 58 59 54;
     54 59 60 55; 
     56 61 62 57;
     57 62 63 58;
     58 63 64 59;
     59 64 65 60];

% Now mirror the structure

nc=node;
nn=size(nc,1); nc=[nc; [nc(:,1) -nc(:,2)]; [-nc(:,1) nc(:,2)]; [-nc(:,1) -nc(:,2)]];
nc=round(1e8*nc)/1e8;
elemnode=[elemnode; elemnode+nn; elemnode+2*nn; elemnode+3*nn];
%There are multiple nodes at intersections (1/4->1) which have to be deleted. 

[~,j]=unique(single(nc),'rows'); node=nc(j,:);
[~,j]=ismember(nc,node,'rows'); %nc-node(j,:)=0
n1=elemnode(:,1); n2=elemnode(:,2); n3=elemnode(:,3); n4=elemnode(:,4); 
nr=[1:size(node,1)]'; nodenr=[1:size(nc,1)]'; m=[nodenr nr(j)];
m1=m(n1,2); m2=m(n2,2); m3=m(n3,2); m4=m(n4,2); clear m
elem =[m1 m2 m3 m4];
if norm( nc(elemnode(:),:)-node(elem(:),:) )>10*eps, error('Error in mesh-function'); end
elemnode=unique(elem,'rows'); nc=node;
nel=size(elemnode,1);
for el=1:nel %To obtain anti-clockwise numbering of nodes
    x=mean(nc(elemnode(el,:),1)); y=mean(nc(elemnode(el,:),2));
    fi=atan2(nc(elemnode(el,:),2)-y,nc(elemnode(el,:),1)-x);
    [~,i]=sort(fi); elemnode(el,:)=elemnode(el,i);
end
 
node=node-min(node);
L=max(node(:,1)); H=max(node(:,2));

% Boundaries
B1=find(node(:,1)<10*eps);
B2=find(node(:,2)>H-10*eps);
B3=find(node(:,1)>L-10*eps);
B4=find(node(:,2)<10*eps);
r=sqrt((0.03-node(:,1)).^2+(0.01-node(:,2)).^2);
BH=find(r<0.00505);

