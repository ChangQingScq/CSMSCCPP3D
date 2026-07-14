function[C] = inverseMapping(curPt, TR1, TR2)
    [T_ID, B] = pointLocation(TR1, curPt);
    curNaN = isnan(T_ID)';
    T_ID = T_ID(~curNaN);
    B = B(~curNaN,:);
    curPt = curPt(~curNaN,:);

    if length(curPt) == 0
        disp('invalid C');
        C=[];
        return;
    end
    %B = cartesianToBarycentric(TR1, T_ID, curPt);
    C = barycentricToCartesian(TR2, T_ID, B);
end