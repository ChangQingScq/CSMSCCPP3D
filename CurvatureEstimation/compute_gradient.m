function [gradientVec, gradientTangent, normalizedNormVec, triangleCenter, A_t] = compute_gradient(vertex, face, scalar)
%%see https://blog.csdn.net/weixin_40211518/article/details/89469816

    node1 = vertex(face(:,1), :);
    node2 = vertex(face(:,2), :);
    node3 = vertex(face(:,3), :);
    A_t = vecnorm(cross(node2- node1, node3- node1, 2),2,2);
    normalizedNormVec = cross(node2- node1, node3- node1, 2)./vecnorm(cross(node2- node1, node3- node1, 2),2,2);
    gradientTangent = (scalar(face(:,2)) - scalar(face(:,1))).* (node1- node3) ./ (2*A_t) + (scalar(face(:,3)) - scalar(face(:,1))).* (node2- node1) ./ (2*A_t);
    gradientVec = cross(normalizedNormVec, gradientTangent, 2);

    triangleCenter = (node1 + node2 + node3) / 3;

end
