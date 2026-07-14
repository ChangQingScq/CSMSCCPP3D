function [energyK] = getEk(v, f, scalar, Bp, Tp, FaceSFMMatrix, boundPtIdx)
%% cal E_Kn


[gradientVec, gradientTangent, normalizedNormVec, triangleCenter, A_t] = compute_gradient(v, f, scalar');
%normGradientVec = gradientVec./vecnorm(gradientVec,2,2);
normGradientTangent = gradientTangent./vecnorm(gradientTangent,2,2);
normGradientVec = gradientVec./vecnorm(gradientVec,2,2);
uComponent = dot(Tp', normGradientTangent')';
vComponent = dot(Bp', normGradientTangent')';
Kn = uComponent.^2.*FaceSFMMatrix(:,1) + 2*uComponent.*vComponent.*FaceSFMMatrix(:,2) + vComponent.^2.*FaceSFMMatrix(:,3);
energyKn = sum((Kn).^2.*A_t);

%% cal E_Kg
FV.faces = f;
FV.vertices = v;
[VertexNormals,Avertex,Acorner,up,vp] = CalcVertexNormals(FV, normalizedNormVec);
Kg = zeros(length(v),1);
for vertexIdx = 1:length(v)
    vOnt = find(any(f == vertexIdx, 2));
    for facesIdx = 1:length(vOnt)
        tNum = vOnt(facesIdx);
        oneFace = circshift(f(tNum,:), 1 - find(vertexIdx==f(tNum,:)));
        a = v(oneFace(1), :); % 当前顶点v1
        b = v(oneFace(2), :); % v2
        c = v(oneFace(3), :); % v3
        bc = v(oneFace(3),:) - v(oneFace(2),:); % edge v3 - v2
        ab = v(oneFace(2),:) - v(oneFace(1),:); % edge v2 - v1
        ac = v(oneFace(3),:) - v(oneFace(1),:); % edge v1 - v3
        cot_beta  = dot(ab, -bc) / (norm(cross(ab, bc)) + eps); % at vertex v2 (中心点)
        cot_gamma = dot(ac, bc) / (norm(cross(ac, bc)) + eps); % at vertex v3
        X = normGradientVec(tNum,:);
        t_Kg = (cot_beta.*dot(ac, X) + cot_gamma.*dot(ab, X));
        Kg(vertexIdx) = Kg(vertexIdx) + t_Kg;
    end
end

%边界处Kg不变,设为0
Kg = - Kg./(2*Avertex);
TR = triangulation(f, v);
boundPtIdx = freeBoundary(TR);
Kg(boundPtIdx) = 0;
%

energyKg = sum((Kg).^2.*Avertex);

%% cal E_K
energyK = energyKn + energyKg;

