function plot_seff(node,elemnode,intpoints,seff)

SE=seff(:); X=intpoints(:,1,:); X=X(:); Y=intpoints(:,2,:); Y=Y(:);

%interpolate value in integration points to nodal points
F = scatteredInterpolant(X,Y,SE,'natural'); se=F(node);

colormap jet
patch('vertices',node,'faces',elemnode,'facevertexcdata',se,'facecolor','interp','linestyle','none');
title('\sigma_e','fontsize',14); axis equal; 
%scaxis([-1 3]); 
axis([min(node(:,1)) max(node(:,1)) min(node(:,2)) max(node(:,2))]); axis off; colorbar; drawnow