function [eta] = closedHarmonicOneFormForSlitMapping(v, f, v_sliced, f_sliced, slice_path)


A = cotLaplacian(v_sliced, f_sliced);


B = [slice_path'; (length(v)+1:length(v_sliced))'];
g = [slice_path'*0+1; (length(v)+1:length(v_sliced))'*0];
nV = size(v_sliced, 1);

I = setdiff(1:nV, B);


AII = A(I,I);
AIB = A(I,B);
ABB = A(B,B);

a = AII \ (I'*0-AIB*g);

gk = zeros(nV,1);
gk(I) = a';
gk(slice_path') = (slice_path)'*0 + 1;

delta_gk = A * gk;
delta_rho = delta_gk(1:length(v));
delta_rho(slice_path) = delta_rho(slice_path) + delta_gk(length(v)+1:end);
% cl2 = mapminmax(rho',0,1)';
% node_colors= cl2.*[0.5 1 0.5]+(1-cl2).*[1 0.5 0];
% trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceVertexCData', node_colors, 'FaceColor', 'interp', 'EdgeColor', 'interp');

%%
A = cotLaplacian(v, f);

rho = A \ -delta_rho;

eta = zeros(length(v_sliced),1);
eta(1:length(rho)) = gk(1:length(rho)) + rho;
eta(length(rho)+1:end) = gk(length(rho)+1:end) + rho(slice_path);


