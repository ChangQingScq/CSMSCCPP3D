% cubic Spline: 
% https://www.cnblogs.com/tiandsp/p/12232392.html

function[re3, rep3, repp3] = splineReparamization(p, n, cornerPt)

%三次均匀b样条
re3=[];
rep3=[];
cornerNum = size(cornerPt, 1);

if 0 == cornerNum
    s = (0 : 1/n : 1-1/n).';
    p = [p; p(1:3,:)];

    for i = 1 : n
        t= s(i) * (length(p)-3) - floor(s(i) * (length(p)-3));

        b0 = 1/6*(1-t)^3;
        b1 = 1/6*(3.*t^3-6*t^2+4);
        b2 = 1/6*(-3*t^3+3*t^2+3*t+1);
        b3 = 1/6*t^3;

        bp0 = -1/2*(1-t)^2;
        bp1 = 3/2*t^2-2*t;
        bp2 = -3/2*t^2+t+1/2;
        bp3 = 1/2*t^2;

        bpp0 = 1-t;
        bpp1 = 3*t-2;
        bpp2 = -3*t+1;
        bpp3 = t;

        idxPt = floor(s(i) * (length(p)-3)) + 1;

        t_re3 = b0*p(idxPt,:)+b1*p(idxPt+1,:)+b2*p(idxPt+2,:)+b3*p(idxPt+3,:);
        t_rep3 = bp0*p(idxPt,:)+bp1*p(idxPt+1,:)+bp2*p(idxPt+2,:)+bp3*p(idxPt+3,:);
        t_repp3 = bpp0*p(idxPt,:)+bpp1*p(idxPt+1,:)+bpp2*p(idxPt+2,:)+bpp3*p(idxPt+3,:);

        re3 = [re3; t_re3];
        rep3 = [rep3; t_rep3*(length(p)-3)/(2*pi)];
        repp3 = [rep3; t_rep3*(length(p)-3)/(2*pi)];
    end
end

if 0 ~= cornerNum
    segmentCurvePt = {};
    for i = 1:cornerNum
        t_segmentCurvePt = [];
        if cornerNum == i
            t_segmentCurvePt = [p(cornerPt(i):end,:); p(1:cornerPt(1),:)];
        else
            t_segmentCurvePt = p(cornerPt(i):cornerPt(i+1),:);
        end
        t_segmentCurvePt = [t_segmentCurvePt(1,:);t_segmentCurvePt(1,:); t_segmentCurvePt; t_segmentCurvePt(end,:); t_segmentCurvePt(end,:)];

        segmentCurvePt{i} = t_segmentCurvePt;
    end

    s = 2*pi*(0 : 1/n*cornerNum : 1-1/n*cornerNum).';
    [s,sp,spp] =  deltw(s,1,3);

    for i = 1 : cornerNum
        for j = 1:length(s)
            t = s(j) * (length(segmentCurvePt{i})-3)/(2*pi) - floor(s(j) * (length(segmentCurvePt{i})-3)/(2*pi));
    
            b0 = 1/6*(1-t)^3;
            b1 = 1/6*(3.*t^3-6*t^2+4);
            b2 = 1/6*(-3*t^3+3*t^2+3*t+1);
            b3 = 1/6*t^3;
    
            bp0 = -1/2*(1-t)^2;
            bp1 = 3/2*t^2-2*t;
            bp2 = -3/2*t^2+t+1/2;
            bp3 = 1/2*t^2;
    
            bpp0 = 1-t;
            bpp1 = 3*t-2;
            bpp2 = -3*t+1;
            bpp3 = t;

            idxPt = floor(s(j) * (length(segmentCurvePt{i})-3)/(2*pi)) + 1;
    
            t_p = segmentCurvePt{i};
            t_re3   = b0  *t_p(idxPt,:) + b1  *t_p(idxPt+1,:) + b2  *t_p(idxPt+2,:) + b3  *t_p(idxPt+3,:);
            t_rep3  = bp0 *t_p(idxPt,:) + bp1 *t_p(idxPt+1,:) + bp2 *t_p(idxPt+2,:) + bp3 *t_p(idxPt+3,:);
            t_repp3 = bpp0*t_p(idxPt,:) + bpp1*t_p(idxPt+1,:) + bpp2*t_p(idxPt+2,:) + bpp3*t_p(idxPt+3,:);

            t_rep3  = t_rep3 * sp(j);
            t_repp3 = t_rep3 * spp(j);

            re3 = [re3; t_re3];

            rep3 = [rep3; t_rep3 * length(segmentCurvePt) * (length(segmentCurvePt{i})-3) / (2*pi)];
            repp3 = [rep3; t_repp3 * length(segmentCurvePt) * (length(segmentCurvePt{i})-3) / (2*pi)];

        end
    end
end
