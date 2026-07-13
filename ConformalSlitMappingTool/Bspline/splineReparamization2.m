function[re3, rep3, repp3] = splineReparamization2(pt, n, cornerPt)

%三次均匀b样条
re3=[];
rep3=[];
repp3=[];
cornerNum = size(cornerPt, 1);

if 0 == cornerNum
    t_query = (0 : 1/n : 1-1/n).' * length(pt);
    [re3, rep3, repp3] = cubic_spline_interp(pt, t_query, "periodic");
    rep3 = rep3 / (length(pt) + 1) * 2*pi;
    repp3 = repp3 / (length(pt) + 1) * 2*pi;
end

if 0 ~= cornerNum
    segmentCurvePt = {};
    for i = 1:cornerNum
        t_segmentCurvePt = [];
        if cornerNum == i
            t_segmentCurvePt = [pt(cornerPt(i):end,:); pt(1:cornerPt(1),:)];
        else
            t_segmentCurvePt = pt(cornerPt(i):cornerPt(i+1),:);
        end
        segmentCurvePt{i} = t_segmentCurvePt;
    end

    s = 2*pi*(0 : 1/n*cornerNum : 1-1/n*cornerNum).';
    [s,sp,spp] =  deltw(s,1,3);

    for i = 1 : cornerNum
        points = segmentCurvePt{i};
        t_query = s / (2*pi) * (length(points) - 1);
        [t_re3, t_rep3, t_repp3] = cubic_spline_interp(points, t_query, "variational");
        t_rep3 = t_rep3 .* sp * (length(points) - 1) / (2*pi);
        t_repp3 = t_repp3 .* spp * (length(points) - 1) / (2*pi);
        re3 = [re3; t_re3];
        rep3 = [rep3; t_rep3];
        repp3 = [repp3; t_repp3];
    end
end
