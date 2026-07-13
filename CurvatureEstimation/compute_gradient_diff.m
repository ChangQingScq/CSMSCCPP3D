function [gradientVecDiff, gradientTangentDiff, normalizedNormVec, triangleCenter, A_t] = compute_gradient_diff(vertex, face, vIdx)

    
    for i = 1:size(face, 1)
        m = find(vIdx == face(i, :));
        face(i, :) = circshift(face(i, :), 2 - m);
    end

    node1 = vertex(face(:,1), :);
    node2 = vertex(face(:,2), :);
    node3 = vertex(face(:,3), :);
    A_t = vecnorm(cross(node2- node1, node3- node1, 2),2,2);
    normalizedNormVec = cross(node2- node1, node3- node1, 2)./vecnorm(cross(node2- node1, node3- node1, 2),2,2);

    gradientTangentDiff = (node1- node3) ./ (2*A_t)  ;
    gradientVecDiff = cross(normalizedNormVec, gradientTangentDiff, 2);

    triangleCenter = (node1 + node2 + node3) / 3;

end