function [re3, rep3, repp3] = cubic_spline_interp(points, t_query, boundMod)
% 三维闭合三次样条插值
% 输入：
%   points : n x 3 的插值点
%   t_query: 要插值的一维参数向量，范围在 [0, n]
%   boundMod: 'periodic'拟合闭合样条曲线 or 'variational'拟合开放样条曲线
% 输出：
%   r  : 插值点坐标
%   r1 : 一阶导数
%   r2 : 二阶导数

n = size(points, 1);
t = 0:n-1;

dim = size(points, 2);
points = [points, zeros(length(points), 3-dim)];

x = points(:, 1);
y = points(:, 2);
z = points(:, 3);
if boundMod == "periodic"
    % 为使曲线闭合，首尾点重复
    x(end+1) = x(1);
    y(end+1) = y(1);
    z(end+1) = z(1);
    t = [t, n];  % 对应的参数向量扩展为 [0, ..., n]
end

% 周期性三次样条拟合
ppx = csape(t, x, boundMod);
ppy = csape(t, y, boundMod);
ppz = csape(t, z, boundMod);

% 计算插值点和导数
re3  = [ppval(ppx, t_query(:)), ppval(ppy, t_query(:)), ppval(ppz, t_query(:))];
rep3 = [ppval(fnder(ppx, 1), t_query(:)), ppval(fnder(ppy, 1), t_query(:)), ppval(fnder(ppz, 1), t_query(:))];
repp3 = [ppval(fnder(ppx, 2), t_query(:)), ppval(fnder(ppy, 2), t_query(:)), ppval(fnder(ppz, 2), t_query(:))];

re3 = re3(:,1:dim);
rep3 = rep3(:,1:dim);
repp3 = repp3(:,1:dim);

end
