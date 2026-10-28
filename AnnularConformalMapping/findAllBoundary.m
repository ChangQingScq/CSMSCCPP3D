function [allBoundary, N] = findAllBoundary(v, f)
TR = triangulation(f, v(:,1), v(:,2), v(:,3));

B = freeBoundary(TR);

allBoundary = {};
startPt = B(1,1);
startPtNum = 1;
for idx = 1:1000
    N = find(B(1:end,2) == startPt);
    allBoundary{idx} = B(startPtNum:N,1);

    if length(B) == N
        break;
    end
    startPtNum = N+1;
    startPt = B(startPtNum,1);
end
