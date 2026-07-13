% author Shen.C.Q 625372005@qq.com 2024.06.30
% modified based on https://www.cnblogs.com/shushen/p/5759679.html

function[fout,vout, freeboundsList, longestBoundIdx] = holeFiller(f, v, model)

[freeboundsList, longestBoundIdx] = getMeshBounds(f, v);

fout = f;
vout = v;
for idx = 1:length(freeboundsList)
    if idx ~= longestBoundIdx
        [fout, vout] = meshFillOneHole(fout, vout, freeboundsList{idx}, model);
    end
end
