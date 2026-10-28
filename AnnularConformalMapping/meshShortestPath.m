function [dist, slice_path, pred] = meshShortestPath(f, v, inner_bd_vid, outer_bd_vid)
% meshShortestPath
% 用 BFS 在三角网格顶点图上求最短路径（无权图）
%
% 输入:
%   f            : nf x 3, 三角面片顶点索引
%   v            : nv x 2 或 nv x 3, 顶点坐标
%   inner_bd_vid : 起点顶点编号
%   outer_bd_vid : 终点顶点编号
%
% 输出:
%   dist         : 最短路径长度（边数）
%   slice_path   : 从 inner_bd_vid 到 outer_bd_vid 的顶点序列
%   pred         : 每个顶点的前驱节点，未访问为 0

    nv = size(v, 1);

    % 构造无向边
    e = [f(:, [1,2]);
         f(:, [2,3]);
         f(:, [3,1])];

    % 邻接矩阵（无权）
    G = sparse(e(:,1), e(:,2), 1, nv, nv);
    G = spones(G + G');   % 转成 0/1 对称邻接矩阵

    [dist, slice_path, pred] = bfsShortestPathSparse(G, inner_bd_vid, outer_bd_vid);
end


function [dist, path, pred] = bfsShortestPathSparse(G, s, t)
% 在稀疏邻接矩阵 G 上，用 BFS 求 s 到 t 的最短路径
%
% 输入:
%   G : nv x nv 稀疏邻接矩阵（无权图，非零表示相邻）
%   s : 起点
%   t : 终点
%
% 输出:
%   dist : 最短距离（边数），不可达时为 Inf
%   path : 最短路径顶点序列，不可达时为空 []
%   pred : 前驱数组

    nv = size(G, 1);

    if s < 1 || s > nv || t < 1 || t > nv
        error('起点或终点索引超出范围。');
    end

    visited = false(nv, 1);
    pred = zeros(nv, 1);
    dist_arr = inf(nv, 1);

    % 简单队列
    queue = zeros(nv, 1);
    head = 1;
    tail = 1;

    queue(tail) = s;
    visited(s) = true;
    dist_arr(s) = 0;

    found = false;

    while head <= tail
        u = queue(head);
        head = head + 1;

        if u == t
            found = true;
            break;
        end

        nbrs = find(G(u, :));
        for k = 1:length(nbrs)
            w = nbrs(k);
            if ~visited(w)
                visited(w) = true;
                pred(w) = u;
                dist_arr(w) = dist_arr(u) + 1;

                tail = tail + 1;
                queue(tail) = w;

                if w == t
                    found = true;
                    break;
                end
            end
        end

        if found
            break;
        end
    end

    if ~visited(t)
        dist = inf;
        path = [];
        return;
    end

    dist = dist_arr(t);

    % 回溯路径
    path = t;
    while path(1) ~= s
        path = [pred(path(1)); path];
    end
end