function[adaptionGrad] = gradientAdaptation(vertexEnergyGradient, dequeIntersecCell, vCSM, freeboundsList, Avertex)

vCSM = [real(vCSM)', imag(vCSM)'];
unitVCSM = vCSM ./ (vecnorm(vCSM, 2, 2) + eps) + eps;
convexIntersect = zeros(length(vCSM), 2);
invConvexIntersect = zeros(length(vCSM), 2);

for vIdx = 1:size(vCSM, 1)
    R0 = vCSM(vIdx, :);
    D = vertexEnergyGradient(vIdx, :);
    P = dequeIntersecCell{vIdx};
    convexIntersect(vIdx, :) = intersectConvexSingleRay(P, R0, D);
    invConvexIntersect(vIdx, :) = intersectConvexSingleRay(P, R0, -D);
end

gradNorm = vecnorm(vertexEnergyGradient, 2, 2);
vecConvexIntersect = convexIntersect - vCSM;
intersectNorm = vecnorm(vecConvexIntersect, 2, 2);

invGradNorm = gradNorm;
vecInvConvexIntersect = invConvexIntersect - vCSM;
invIntersectNorm = vecnorm(vecInvConvexIntersect, 2, 2);

boundIterVec = zeros(length(freeboundsList), 1);
ptOnBound = [];
for bdIdx = 1:length(freeboundsList)
    limitDistPt = freeboundsList{bdIdx};
    boundIterVec(bdIdx) = sum(dot(vertexEnergyGradient(limitDistPt, :), ...
                          unitVCSM(limitDistPt, :), 2) .* Avertex(limitDistPt, :) ) ./ ...
                          sum(Avertex(limitDistPt, :));
    ptOnBound = [ptOnBound; freeboundsList{bdIdx}];
end
ptNotOnBound = setdiff((1:length(vCSM))',ptOnBound);

%%
adaptionGrad = vecConvexIntersect;

if 1

normRatio = gradNorm(ptNotOnBound)./(intersectNorm(ptNotOnBound) + eps);
maxNormRatio = max(normRatio);

quantile_persent_pt = ceil(1/3 * length(normRatio));
quantile_pt = sort(normRatio);
normRatioQ_1_3 = quantile_pt(quantile_persent_pt); %normRatioQ_1_3 normRatio的1/3位数
lambda = normRatioQ_1_3/(-log(1 - 1/5));
adaptationRatio = 1 - exp(-normRatio/(lambda + eps));

avgNormRatio = sum(normRatio .* Avertex(ptNotOnBound) ./ sum(Avertex(ptNotOnBound)));
adaptationRatio = normRatio ./ (normRatio + avgNormRatio + eps);
adaptionGrad(ptNotOnBound, :) = vertexEnergyGradient(ptNotOnBound, :) ./(normRatio + eps) .* adaptationRatio;
adaptionGrad(ptNotOnBound, :) = vertexEnergyGradient(ptNotOnBound, :) ./(normRatio + eps);
%adaptionGrad(ptNotOnBound, :) = vertexEnergyGradient(ptNotOnBound, :) ./maxNormRatio;
% invNormRatio = invGradNorm(ptNotOnBound)./(invIntersectNorm(ptNotOnBound) + eps);
% maxInvNormRatio = max(invNormRatio);
% 
% expSigma = -1/ maxInvNormRatio * log(1/reachPersent - 1);
% adaptationInvRatio = 1 ./ (1 + exp(-expSigma * invNormRatio));
% adaptionInvGrad(ptNotOnBound) = vecInvConvexIntersect(ptNotOnBound) .* adaptationInvRatio;
end


% adaptionGrad(ptNotOnBound, :) = vertexEnergyGradient(ptNotOnBound, :) / maxRatio;
% adaptionGrad(ptNotOnBound, :) = vertexEnergyGradient(ptNotOnBound, :) / maxRatio;
% adaptionGrad(ptNotOnBound, :) = vertexEnergyGradient(ptNotOnBound, :)./ gradNorm(ptNotOnBound, :)...
%                                 .* intersectNorm(ptNotOnBound, :);

if 0
reachPersent = 1; % 边界统一沿经向的运动达到允许最大运动距离的20%
for bdIdx = 1:length(freeboundsList)
%     if bdIdx == 2 || bdIdx == 3
%         reachPersent = 1;
%     else
%         reachPersent = 1;
%     end
    limitDistPt = freeboundsList{bdIdx};

    radialVecNorm = vecnorm(vecConvexIntersect(limitDistPt, :) + vCSM(limitDistPt, :),2,2);
    dotVar = dot(vecConvexIntersect(limitDistPt, :),unitVCSM(limitDistPt, :),2);
    crossOrigin = find((dotVar + vecnorm(vCSM(limitDistPt, :), 2, 2)) < 0);
    radialVecNorm(crossOrigin) = vecnorm(vCSM(limitDistPt(1), :), 2, 2) * 0.5;
    radialVec = vCSM(limitDistPt, :) ./ (vecnorm(vCSM(limitDistPt, :), 2, 2) + eps) .* radialVecNorm;

    invRadialVecNorm = vecnorm(vecInvConvexIntersect(limitDistPt, :) + vCSM(limitDistPt, :),2,2);
    invDotVar = dot(vecInvConvexIntersect(limitDistPt, :),unitVCSM(limitDistPt, :),2);
    crossOrigin = find((invDotVar + vecnorm(vCSM(limitDistPt, :), 2, 2)) < 0);
    invRadialVecNorm(crossOrigin) = vecnorm(vCSM(limitDistPt(1), :), 2, 2) * 0.5;
    invRadialVec = vCSM(limitDistPt, :) ./ (vecnorm(vCSM(limitDistPt, :), 2, 2) + eps) .* invRadialVecNorm;
 
    sortRadial = sort([dotVar, invDotVar], 2);
    maxDiff = 3; %正反向之差不超过maxDiff倍
    rMove = 0;
    r = 0;
    startPt = [vCSM(limitDistPt, :), vecConvexIntersect(limitDistPt, :)];
    if boundIterVec(bdIdx) > 0
        rMove = min([sortRadial(:,2); maxDiff * min(abs(sortRadial(:,1)))]);
        r = mean(vecnorm(vCSM(limitDistPt, :),2,2)) + rMove;
        startPt(find(dotVar > 0), 3:4) = radialVecNorm(find(dotVar > 0), :) - vCSM(limitDistPt(find(dotVar > 0)), :);
        startPt(find(dotVar <= 0), 3:4) = invRadialVecNorm(find(dotVar <= 0), :) - vCSM(limitDistPt(find(dotVar <= 0)), :);
    else
        rMove = max([sortRadial(:,1); -maxDiff * min(abs(sortRadial(:,2)))]);
        r = mean(vecnorm(vCSM(limitDistPt, :),2,2)) + rMove;
        startPt(find(dotVar > 0), 3:4) = invRadialVecNorm(find(dotVar > 0), :) - vCSM(limitDistPt(find(dotVar > 0)), :);
        startPt(find(dotVar <= 0), 3:4) = radialVecNorm(find(dotVar <= 0), :) - vCSM(limitDistPt(find(dotVar <= 0)), :);
    end
    [P1, ~] = intersectLineCircle(r, rMove, startPt);
    adaptionGrad(limitDistPt, :) = P1 - vCSM(limitDistPt, :);
    if 0 == rMove
        adaptionGrad(limitDistPt, :) = 0;
    end
    


% 
%     normVecConvexIntersectOnBd = vecConvexIntersect(limitDistPt, :)./vecnorm(vecConvexIntersect(limitDistPt, :), 2, 2);
%     adaptionGrad(limitDistPt, :) = vecConvexIntersect(limitDistPt, :) ./ dotVar * rMove;
%     
%     adaptionGrad(limitDistPt, :) * (rMove + 1);
% 
% 
% 
%     sortRadial = sort([radialVec, invRadialVec], 2);
%     curR = mean(vCSM(limitDistPt, :));
% 
%     r = 0;
%     startPt = [vCSM(limitDistPt, :), vecConvexIntersect(limitDistPt, :)];
%     if boundIterVec(bdIdx) > 0
%         r = min(sortRadial(:, 2));
%         if r>max(abs(vCSM(limitDistPt, :)))*1.1
%             r=max(vecnorm(vCSM(limitDistPt, :),2,2))*1.1;
%         end
%         
%         startPt(find(dotVar > 0), 3:4) = vecConvexIntersect(limitDistPt(find(dotVar > 0)), :);
%         startPt(find(dotVar <= 0), 3:4) = -vecConvexIntersect(limitDistPt(find(dotVar <= 0)), :);
%     else
%         r = max(sortRadial(:, 1));
%         if r>max(abs(vCSM(limitDistPt, :)))*1.1
%             r=max(vecnorm(vCSM(limitDistPt, :),2,2))*1.1;
%         end
%         startPt(find(dotVar > 0), 3:4) = -vecConvexIntersect(limitDistPt(find(dotVar > 0)), :);
%         startPt(find(dotVar <= 0), 3:4) = vecConvexIntersect(limitDistPt(find(dotVar <= 0)), :);
%     end
% 
%     [P1, ~] = intersectLineCircle(r, startPt);
%     adaptionGrad(limitDistPt, :) = P1 - vCSM(limitDistPt, :);

end
end

reachPersent = 10;
for bdIdx = 1:length(freeboundsList)

    limitDistPt = freeboundsList{bdIdx};

    radialVec = vecConvexIntersect(limitDistPt, :) + vCSM(limitDistPt, :);
    dotVar = dot(vecConvexIntersect(limitDistPt, :), unitVCSM(limitDistPt, :), 2);
    crossOrigin = find((dotVar + vecnorm(vCSM(limitDistPt, :), 2, 2)) < 0.5 * vecnorm(vCSM(limitDistPt(1), :), 2, 2));
    radialVec(crossOrigin, :) = vCSM(limitDistPt(crossOrigin), :) * 0.5;
    dotVar(crossOrigin) = -0.5 * vecnorm(vCSM(limitDistPt(crossOrigin), :), 2, 2);
    crossOrigin = find((dotVar + vecnorm(vCSM(limitDistPt, :), 2, 2)) > 1.1 * vecnorm(vCSM(limitDistPt(1), :), 2, 2));
    radialVec(crossOrigin, :) = vCSM(limitDistPt(crossOrigin), :) * 1.1;
    dotVar(crossOrigin) = 0.1 * vecnorm(vCSM(limitDistPt(crossOrigin), :), 2, 2);

%     radialVec = vecConvexIntersect(limitDistPt, :) + vCSM(limitDistPt, :);
%     dotVar = dot(vecConvexIntersect(limitDistPt, :),unitVCSM(limitDistPt, :),2);
%     crossOrigin = find((dotVar + vecnorm(vCSM(limitDistPt, :), 2, 2)) < 0);
%     vecConvexIntersectUnit = vecConvexIntersect(limitDistPt(crossOrigin), :) ./ vecnorm(vecConvexIntersect(limitDistPt(crossOrigin), :), 2, 2);
%     radialVec(crossOrigin) = vecConvexIntersectUnit * 0.5 * vecnorm(vCSM(limitDistPt(crossOrigin), :), 2, 2) + vCSM(limitDistPt(crossOrigin), :);
%     dotVar(crossOrigin) = -0.5 * vecnorm(vCSM(limitDistPt(crossOrigin), :), 2, 2);

    invRadialVec = vecInvConvexIntersect(limitDistPt, :) + vCSM(limitDistPt, :);
    invDotVar = dot(vecInvConvexIntersect(limitDistPt, :), unitVCSM(limitDistPt, :), 2);
    crossOrigin = find((invDotVar + vecnorm(vCSM(limitDistPt, :), 2, 2)) < 0.5 * vecnorm(vCSM(limitDistPt(1), :), 2, 2));
    invRadialVec(crossOrigin, :) = vCSM(limitDistPt(crossOrigin), :) * 0.5;
    invDotVar(crossOrigin) = -0.5 * vecnorm(vCSM(limitDistPt(crossOrigin), :), 2, 2);
    crossOrigin = find((invDotVar + vecnorm(vCSM(limitDistPt, :), 2, 2)) > 1.1 * vecnorm(vCSM(limitDistPt(1), :), 2, 2));
    invRadialVec(crossOrigin, :) = vCSM(limitDistPt(crossOrigin), :) * 1.1;
    invDotVar(crossOrigin) = 0.1 * vecnorm(vCSM(limitDistPt(crossOrigin), :), 2, 2);


        %dr = vecnorm([radialVec; invRadialVec], 2, 2) - mean(vecnorm(vCSM(limitDistPt, :), 2, 2));
    
    sortRadial = sort([dotVar, invDotVar], 2);
    maxDiff = 3; %正反向之差不超过maxDiff倍
    rMove = 0;
    r = 0;
    startPt = [vCSM(limitDistPt, :), vecConvexIntersect(limitDistPt, :)];
    if boundIterVec(bdIdx) > 0
        rMove = min([sortRadial(:,2); maxDiff * min(abs(sortRadial(:,1)))]);
        r = mean(vecnorm(vCSM(limitDistPt, :),2,2)) + rMove;
        startPt(find(dotVar > 0), 3:4) = radialVec(find(dotVar > 0), :) - vCSM(limitDistPt(find(dotVar > 0)), :);
        startPt(find(dotVar <= 0), 3:4) = invRadialVec(find(dotVar <= 0), :) - vCSM(limitDistPt(find(dotVar <= 0)), :);
    else
        rMove = max([sortRadial(:,1); -maxDiff * min(abs(sortRadial(:,2)))]);
        r = mean(vecnorm(vCSM(limitDistPt, :),2,2)) + rMove;
        startPt(find(dotVar > 0), 3:4) = invRadialVec(find(dotVar > 0), :) - vCSM(limitDistPt(find(dotVar > 0)), :);
        startPt(find(dotVar <= 0), 3:4) = radialVec(find(dotVar <= 0), :) - vCSM(limitDistPt(find(dotVar <= 0)), :);
    end
    [P1, ~] = intersectLineCircle(r, rMove, startPt);
    adaptionGrad(limitDistPt, :) = P1 - vCSM(limitDistPt, :);

%     if bdIdx == 2 || bdIdx == 3
%         reachPersent = 20;
%     else
%         reachPersent = 5;
%     end
    adaptionGrad(limitDistPt, :) = adaptionGrad(limitDistPt, :) * reachPersent;
    if 0 == rMove
        adaptionGrad(limitDistPt, :) = 0;
    end
end

