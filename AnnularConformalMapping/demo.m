clear
close all
clc
[v, f] = generatingMeshForAnnularMapping();

[v_sliced, f_sliced, f0, eta, lambda] = annularConformalMapping(v, f);
%lambda= 1/lambda;%lambda=1;
%lambda = lambda+1;
TR = triangulation(f_sliced, f0, eta*lambda);

figure()
trisurf(f_sliced, f0, eta*lambda, f0*0, 'FaceVertexCData', eta*lambda, 'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal

color = round((eta*lambda)/0.11444)+round((f0)/0.11444);
color = mod(color,2);
color = color*[1 1 1];

figure()
trisurf(f_sliced, f0, eta*lambda, f0*0, 'FaceVertexCData', color, 'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal

figure()
trisurf(f_sliced, v_sliced(:,1), v_sliced(:,2), v_sliced(:,3), 'FaceVertexCData', color, 'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal
hold on
plot3(v_sliced(length(v)+1:end,1), v_sliced(length(v)+1:end,2), v_sliced(length(v)+1:end,3))

