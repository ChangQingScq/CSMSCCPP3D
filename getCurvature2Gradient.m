function [curvature2Gradient, gradientVec, triangleCenter, A_t] = getCurvature2Gradient(v, f, scalar, Bp, Tp, FaceSFMMatrix)

[gradientVec, gradientTangent, normalizedNormVec, triangleCenter, A_t] = compute_gradient(v, f, scalar);

normGradientVec = gradientVec./vecnorm(gradientVec,2,2);

uComponent = dot(Tp', normGradientVec')';
vComponent = dot(Bp', normGradientVec')';

%curvature2Gradient = zeros(length(FaceSFMMatrix),1);

curvature2Gradient = uComponent.^2.*FaceSFMMatrix(:,1) + 2*uComponent.*vComponent.*FaceSFMMatrix(:,2) + vComponent.^2.*FaceSFMMatrix(:,3);

% for idx = 1:length(FaceSFM)
%     curFaceSFM = cell2mat(FaceSFM(idx));
%     currcurvature2Gradient = uComponent(idx).^2.*curFaceSFM(1,1) + 2*uComponent(idx).*vComponent(idx).*curFaceSFM(1,2) + vComponent(idx).^2.*curFaceSFM(2,2);
%     curvature2Gradient(idx) = currcurvature2Gradient;
% end

%[EFG,u_vec,v_vec] = compute_second_basic_form(v, f);


