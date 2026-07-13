function cornerIdx = cornerIdentify(v, freeboundsList)
    % to do
    cornerIdx = {};
    for idx = length(freeboundsList)
        cornerIdx{idx} = [];
        if idx == 5
            %cornerIdx{idx} = [6 9 19 24]';
            cornerIdx{idx} = extractSharpCorners3D(v(freeboundsList{idx},:),4);
        end
    end
end

function minAnglePoints = extractSharpCorners3D(points, numCorners)
    % 输入：points - 封闭曲线的离散点，N x 3矩阵
    % 输入：numCorners - 提取的最小角度点的数量
    
    % 点的数量
    n = size(points, 1);
    
    % 计算每个点与相邻两点构成的角度
    angles = zeros(n, 1);  % 存储每个点的角度

    prev_point = circshift(points,  1);    % 前一个点
    next_point = circshift(points, -1);    % 后一个点


    % 计算向量
    v1 = prev_point - points;  % 前向量
    v2 = next_point - points;  % 后向量
    % 计算叉积，得到平面法向量
    cross_prod = cross(v1, v2);
    norm_cross = vecnorm(cross_prod,2,2);
    angles = acos(dot(v1', v2')' ./ (vecnorm(v1,2,2) .* vecnorm(v2,2,2)));  % 计算夹角
    % 找到角度最小的若干个点
    [~, sortedIndices] = sort(angles);
    minAnglePoints = sortedIndices(1:numCorners);  % 提取角度最小的点
    minAnglePoints = sort(minAnglePoints);

end

