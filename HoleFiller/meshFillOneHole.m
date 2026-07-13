% author Shen.C.Q 625372005@qq.com 2024.06.30
% modified based on:
% https://www.cnblogs.com/shushen/p/5759679.html
% [1] Olga Sorkine and Daniel Cohen-Or. 2004. Least-Squares Meshes. In Proceedings of the Shape Modeling International 2004 (SMI '04). IEEE Computer Society, Washington, DC, USA, 191-199.

function [fout, vout] = meshFillOneHole(f, v, bound, model)

fout = f;
vout = v;

if model == "NotFill"
    return;
end

if model == "BaryCenter"
    baryCenter = sum(v(bound,:))/length(bound);
    vout = [vout; baryCenter];
    new_f = [bound, circshift(bound,1), (bound*0+ length(v) +1)];
    fout = [fout; new_f];
    return;
end

%T = triangulation(f, v);
% hold on
% axis equal
% trimesh(T);
% myPlot(v(bound,:),'g-*')
t_bound = bound;
avgfreeboundsLength = 0;
whileIDx = 0;
while 1
    whileIDx = whileIDx + 1;
    if 3 == length(bound)
        fout = [fout; bound'];
        bound = [];
        break;
    end

    boundShift = [bound(2:end);bound(1)];
    boundPt = vout(bound,:);
    boundShiftPt = vout(boundShift, :);
    if 1 == whileIDx
        avgfreeboundsLength = sum(vecnorm(boundPt - boundShiftPt, 2, 2))./length(boundPt);
    end
    boundVec = boundShiftPt - boundPt;
    boundVecShift = -[ boundVec(end, :); boundVec(1:end-1,:)];

    cos_angle = dot(boundVec', boundVecShift')' ./ (vecnorm(boundVec, 2, 2) .* vecnorm(boundVecShift, 2, 2));
    angle_rad = acos(cos_angle);
    patchFillIdx = find(angle_rad == min(angle_rad),1);
    patchFillIdxPrev = 0;
    patchFillIdxNext = 0;
    if patchFillIdx == 1
        patchFillIdxPrev = length(bound);
    else
        patchFillIdxPrev = patchFillIdx-1;
    end
    if patchFillIdx == length(bound)
        patchFillIdxNext = 1;
    else
        patchFillIdxNext = patchFillIdx+1;
    end

    edgeLength = norm(vout(bound(patchFillIdxPrev),:) - vout(bound(patchFillIdxNext),:));
    if edgeLength < avgfreeboundsLength/2
        add_f = [bound(patchFillIdxPrev), bound(patchFillIdx), bound(patchFillIdxNext)];
        fout = [fout; add_f];

        addEdge = vout( bound([patchFillIdxPrev,patchFillIdxNext],:),:);
        %plot3(addEdge(:,1),addEdge(:,2),addEdge(:,3),'r-*')
        bound(patchFillIdx) = [];
    else
        add_v = (vout(bound(patchFillIdxPrev),:) + vout(bound(patchFillIdxNext),:))/2;
        vout = [vout; add_v];
        add_f1 = [bound(patchFillIdxPrev), length(vout), bound(patchFillIdx)];
        add_f2 = [length(vout), bound(patchFillIdx), bound(patchFillIdxNext)];
        fout = [fout; add_f1; add_f2];
        bound(patchFillIdx) = length(vout);
        addEdge = vout( [bound(patchFillIdxPrev),end,bound(patchFillIdxNext)],: );
        %plot3(addEdge(:,1),addEdge(:,2),addEdge(:,3),'b-*')
    end
end



if model == "LSM"

    v_connected = {};
    fillMesh = [length(v)+1:length(vout),t_bound']';
    for idx = 1:length(fillMesh)
        v_connected_num = [];
        for idx2 = 1:3
            v_connected_num = [v_connected_num; find(fillMesh(idx) == fout(length(f)+1:end, idx2))];
        end
        v_connected_num = fout(v_connected_num + length(f), :);
        v_connected{idx} = unique(v_connected_num);
    end

    %L = zeros(length(fillMesh));
    L = sparse(length(fillMesh));
    for idx = 1:length(fillMesh)
        connectedPt = v_connected{idx};
        d = length(connectedPt)-1;
        for idx2 = 1:length(connectedPt)
            connectedPtNum = find(fillMesh == connectedPt(idx2));
            L(idx, connectedPtNum) = -1/d;
        end
        L(idx, idx) = 1;
    end
    
    %P = zeros(length(t_bound), length(fillMesh));
    P = sparse(length(t_bound), length(fillMesh));
    for idx = 1:length(t_bound)
        connectedPt = t_bound(idx);
        connectedPtNum = find(fillMesh == connectedPt);
        P(idx, connectedPtNum) = 1;
    end
    
    A = [L; P];
    
    B = [zeros(length(fillMesh),3); vout(t_bound, :)];
    
    %vLSM = inv(A'*A)*A'*B;
    vLSM = (A'*A)\A'*B;

    vout(fillMesh, :) = vLSM;
end

if model == "RBF"
% to do
end


