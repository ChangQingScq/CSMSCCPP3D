function [W] = edgeProductWeight(V, F)
% author Qing Mu 625372005@qq.com
% 

nV = size(V, 1);

allEdge = [];
for idx = 1:length(F)
    allEdge = [allEdge; F(idx,1), F(idx,2); F(idx,2), F(idx,3); F(idx,3), F(idx,1)];
end
allEdge = sort(allEdge,2);
allEdge = unique(allEdge,"rows");

allTOnE = zeros(length(allEdge), 8);
for idx = 1:length(allEdge)
    [c1,r1] = find(F == allEdge(idx, 1));
    [c2,r2] = find(F == allEdge(idx, 2));
    T_on_E = intersect(c1, c2);

    if 1 == length(T_on_E)
        allTOnE(idx,1:3) = F(T_on_E(1),:);
        allTOnE(idx,4) = setdiff(F(T_on_E(1),:), allEdge(idx, 1:2));
    elseif 2 == length(T_on_E)
        allTOnE(idx,1:3) = F(T_on_E(1),:);
        allTOnE(idx,4) = setdiff(F(T_on_E(1),:), allEdge(idx, 1:2));
        allTOnE(idx,5:7) = F(T_on_E(2),:);
        allTOnE(idx,8) = setdiff(F(T_on_E(2),:), allEdge(idx, 1:2));
    end
end

W = sparse(nV,nV);

edge_cot = zeros(length(allEdge),2);
for idx = 1:length(allTOnE)
    edge1 = V(allEdge(idx, 1) ,:) - V(allTOnE(idx, 4) ,:);
    edge2 = V(allEdge(idx, 2) ,:) - V(allTOnE(idx, 4) ,:);
    Cos = dot(edge1, edge2)/vecnorm(edge1)/vecnorm(edge2);
    cot = cotd(rad2deg(acos(Cos)));
    edge_cot(idx, 1) = cot;
    if 0 ~= allTOnE(idx, 8)
        edge3 = V(allEdge(idx, 2) ,:) - V(allTOnE(idx, 8) ,:);
        edge4 = V(allEdge(idx, 1) ,:) - V(allTOnE(idx, 8) ,:);
        Cos = dot(edge3, edge4)/vecnorm(edge3)/vecnorm(edge4);
        cot = cotd(rad2deg(acos(Cos)));
        edge_cot(idx, 2) = cot;
    end
    W(allEdge(idx,1), allEdge(idx,2)) = edge_cot(idx, 1) + edge_cot(idx, 2);
    W(allEdge(idx,2), allEdge(idx,1)) = edge_cot(idx, 1) + edge_cot(idx, 2);
end



