function [sx,sy,sxy]=interpolate2nodes(s,intpoints,node)
% [Sx,Sy,Sxy]=interpolate2nodes(s,intpoints,node)
%
% function to interpolate variable s (e.g. stresses or strains)
% in element's integration points to global system nodes
%
% Example:
% [sx,sy,sxy]=interpolate2nodes(stress,intpoints,node);
%
% [sx,sy,sxy] is the stress components at nodal positions
%
% Per Isaksson 2020-10-07

%Extract values from multi-dimensional array
SX =s(:,1,:); SX =SX(:);
SY =s(:,2,:); SY =SY(:);
SXY=s(:,3,:); SXY=SXY(:);
X=intpoints(:,1,:); X=X(:);
Y=intpoints(:,2,:); Y=Y(:);

F = scatteredInterpolant(X,Y,SX,'natural'); sx=F(node);
F = scatteredInterpolant(X,Y,SY,'natural'); sy=F(node);
F = scatteredInterpolant(X,Y,SXY,'natural'); sxy=F(node);
