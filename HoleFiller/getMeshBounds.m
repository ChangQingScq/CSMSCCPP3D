function[freeboundsList, longestBoundIdx] = getMeshBounds(f, v)
T = triangulation(f, v);
freebounds = freeBoundary(T); %initialBounds

freeboundsIdx = [];
for idx = 1:length(freebounds)-1
    if freebounds(idx, 2) ~= freebounds(idx + 1, 1)
        freeboundsIdx = [freeboundsIdx; idx];
    end
end

if length(freeboundsIdx) == 0
    freeboundsList{1} = freebounds(:, 1);
    longestBoundIdx = 1;
    return;
end

freeboundsList = {};
for idx = 1:length(freeboundsIdx)+1
    if 1 == idx
        freeboundsList{idx} = freebounds(1:freeboundsIdx(idx), 1);
    elseif length(freeboundsIdx) +1 == idx
        freeboundsList{idx} = freebounds(freeboundsIdx(idx-1)+1 :end, 1);
    else
        freeboundsList{idx} = freebounds(freeboundsIdx(idx-1) + 1 : freeboundsIdx(idx), 1);
    end
end

freeboundsLengthList = zeros(length(freeboundsList),1);
for idx = 1:length(freeboundsList)
    bound = freeboundsList{idx};
    boundShift = [bound(2:end);bound(1)];
    freeboundsLengthList(idx) = sum(vecnorm(v(bound,:) - v(boundShift,:), 2, 2));    
end
longestBoundIdx = find(freeboundsLengthList == max(freeboundsLengthList));
