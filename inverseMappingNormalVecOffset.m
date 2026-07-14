function[C, triNorm] = inverseMappingNormalVecOffset(curPt, TR1, TR2)
    FN = vertexNormal(TR2); %顶点法矢
    T_ID = pointLocation(TR1, curPt);
    curNaN = isnan(T_ID)';
    T_ID = T_ID(~curNaN);
    curPt = curPt(~curNaN,:);
    if length(curPt) == 0
        disp('invalid C');
        C=[];
        return;
    end
    B = cartesianToBarycentric(TR1, T_ID, curPt);
    C = barycentricToCartesian(TR2, T_ID, B);

    triNode = TR2.ConnectivityList(T_ID,:);
    triNorm = FN(triNode(:,1),:).*B(:,1) + FN(triNode(:,2),:) .*B(:,2) + FN(triNode(:,3),:) .*B(:,3);
    triNorm = triNorm./vecnorm(triNorm, 2, 2);
    triNorm = -triNorm;
end