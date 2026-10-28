function [eta, v_sliced, f_sliced, slice_path] = closedHarmonicOneForm(v, f)


bd = meshboundaries(f); % get all boundaries
nb = length(bd); % number of boundaries
if nb ~= 2
    error('The input mesh is not topologically equivalent to annulus!');
end

% determine which boundary is the outer one
bd1 = bd{1};
bd2 = bd{2};
if size(v,2) == 2
    % use the Pick's formula
    area1 = abs(sum(v(bd1,1).*v(bd1([2:length(bd1),1]),2)) - ...
        sum(v(bd1([2:length(bd1),1]),1).*v(bd1,2)));
    area2 = abs(sum(v(bd2,1).*v(bd2([2:length(bd2),1]),2)) - ...
        sum(v(bd2([2:length(bd2),1]),1).*v(bd2,2)));
    if area1 > area2
        outer_bd_id = 1;
        inner_bd_id = 2;
    else
        outer_bd_id = 2;
        inner_bd_id = 1;
    end        
else
    % pick the longer boundary as the outer boundary
    bd1_length = sum(sqrt(sum((v(bd1,:) - v(bd1([2:end,1]),:)).^2,2)));
    bd2_length = sum(sqrt(sum((v(bd2,:) - v(bd2([2:end,1]),:)).^2,2)));
    if bd1_length > bd2_length
        outer_bd_id = 1;
        inner_bd_id = 2;
    else 
        outer_bd_id = 2;
        inner_bd_id = 1;
    end
end
%% slice the mesh 
e = [f(:,[1,2]);f(:,[2,3]);f(:,[3 1])];
nv = length(v);
ne = length(e);

% arbitrarily choose a vertex at the inner boundary
inner_bd_vid = bd{inner_bd_id}(1);

% find the closest vertex at the outer boundary
[~,outer_bd_vid_temp] = min(sum((v(bd{outer_bd_id},:) - repmat(v(inner_bd_vid,:),length(bd{outer_bd_id}),1)).^2,2));
outer_bd_vid = bd{outer_bd_id}(outer_bd_vid_temp);
% find the shortest path between the two vertices
G = sparse(e(:,1),e(:,2),ones(ne,1),nv,nv);
G = (G + G');
[~, slice_path, ~] = graphshortestpath(G,inner_bd_vid,outer_bd_vid);

% ensure no other vertices in slice_path are on the boundaries
[~,id1,~] = intersect(slice_path,bd{inner_bd_id});
[~,id2,~] = intersect(slice_path,bd{outer_bd_id});
slice_path = slice_path(max(id1):min(id2));

% slice mesh along the shortest path 
[v_sliced,f_sliced] = slice_mesh_annulus(v,f,slice_path,bd{inner_bd_id});

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

A = cotLaplacian(v, f);

rho = A \ -delta_rho;

eta = zeros(length(v_sliced),1);
eta(1:length(rho)) = gk(1:length(rho)) + rho;
eta(length(rho)+1:end) = gk(length(rho)+1:end) + rho(slice_path);


