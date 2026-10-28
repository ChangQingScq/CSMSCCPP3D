function[all_v_sliced, all_f_sliced, all_slice_path, outer_bd_id] = sliceMesh(v, f)
% modified by Shen Changqing
%
% Slice multiply connected mesh into annulus patches
%
% The shortest slit path is computed on a reduced mesh:
% triangles connected with other boundaries are removed.
%
% The original mesh is used for final slicing.

%% =========================================================
% Find boundaries
%% =========================================================

bd = meshboundaries(f);
nb = length(bd);


%% =========================================================
% Find outer boundary
%% =========================================================

outer_bd_id = 0;
longest_bd_length = 0;

for idx = 1:length(bd)

    t_bd = bd{idx};
    t_bd_shift = circshift(t_bd,-1);

    t_bd_length = sum(vecnorm(...
        v(t_bd,:) - v(t_bd_shift,:),2,2));

    if t_bd_length > longest_bd_length
        longest_bd_length = t_bd_length;
        outer_bd_id = idx;
    end

end


%% =========================================================
% Slice mesh
%% =========================================================

all_v_sliced = {};
all_f_sliced = {};
all_slice_path = {};


for inner_bd_id = 1:length(bd)


    % outer boundary does not need slicing
    if inner_bd_id == outer_bd_id

        all_v_sliced{inner_bd_id}=0;
        all_f_sliced{inner_bd_id}=0;
        all_slice_path{inner_bd_id}=0;

        continue;

    end



    %% =====================================================
    % Remove other boundary connected triangles
    %% =====================================================

    other_bd_id = setdiff(1:length(bd),...
                    [outer_bd_id,inner_bd_id]);


    other_bd_vertices = [];

    for k = other_bd_id
        other_bd_vertices = [
            other_bd_vertices;
            bd{k}(:)
            ];
    end

    other_bd_vertices = unique(other_bd_vertices);



    % remove faces containing other boundary vertices
    remove_face = any(ismember(f,...
                    other_bd_vertices),2);


    f_graph = f(~remove_face,:);



    %% =====================================================
    % Construct graph on reduced mesh
    %% =====================================================

    e = [
        f_graph(:,[1,2]);
        f_graph(:,[2,3]);
        f_graph(:,[3,1])
        ];


    % remove duplicate edges
    e = sort(e,2);
    e = unique(e,'rows');


    % edge length as weight
    edge_length = vecnorm(...
        v(e(:,1),:) - v(e(:,2),:),2,2);


    nv = size(v,1);


    G = sparse(...
        e(:,1),...
        e(:,2),...
        edge_length,...
        nv,...
        nv);


    G = G + G';



    %% =====================================================
    % Find closest vertex pair between two boundaries
    %% =====================================================

    inner_vertices = bd{inner_bd_id};
    outer_vertices = bd{outer_bd_id};


    dist_mat = pdist2(...
        v(inner_vertices,:),...
        v(outer_vertices,:));


    [~,min_id] = min(dist_mat(:));


    [inner_idx,outer_idx] = ind2sub(...
        size(dist_mat),min_id);


    inner_bd_vid = inner_vertices(inner_idx);
    outer_bd_vid = outer_vertices(outer_idx);



    %% =====================================================
    % Shortest path
    %% =====================================================

    [dist_path,slice_path] = ...
        graphshortestpath(...
        G,...
        inner_bd_vid,...
        outer_bd_vid);



    % no valid path
    if isempty(slice_path) || isinf(dist_path)

        warning(['No valid slit path between boundary ',...
            num2str(inner_bd_id)]);

        all_v_sliced{inner_bd_id}=0;
        all_f_sliced{inner_bd_id}=0;
        all_slice_path{inner_bd_id}=0;

        continue;

    end



    %% =====================================================
    % Slice original mesh
    %% =====================================================

    [v_sliced,f_sliced] = ...
        slice_mesh_annulus(...
        v,...
        f,...
        slice_path,...
        bd{inner_bd_id});


    all_v_sliced{inner_bd_id}=v_sliced;
    all_f_sliced{inner_bd_id}=f_sliced;
    all_slice_path{inner_bd_id}=slice_path;


end

end