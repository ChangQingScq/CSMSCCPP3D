function[basicEnergy] = calBasicEnergy(gradientVec, k_c, curvature2Gradient, A_t, EsMode)

normGradientVec = (vecnorm(gradientVec,2,2) + eps);
uniformCharacter = ((k_c + curvature2Gradient)/8 + eps);

switch (EsMode)
    case 1 % zou qiang 2014
        basicEnergy = sum ((normGradientVec-sqrt(uniformCharacter)).^2 .*A_t);
    case 2 % MIPS
        basicEnergy = sum ( ...
            (((k_c + curvature2Gradient)/8).^0.5./normGradientVec + normGradientVec./(k_c + curvature2Gradient)/8).^0.5.*A_t);
    case 3
        basicEnergy = sum ( ...
            (((k_c + curvature2Gradient)/8)./normGradientVec.^2 + normGradientVec.^2./(k_c + curvature2Gradient)/8).*A_t);
    case 4
        Average = sum ( ...
            uniformCharacter./normGradientVec.^2 .*A_t)./sum(A_t);
        basicEnergy = sum ( ...
            (uniformCharacter./normGradientVec./Average + Average.*normGradientVec./uniformCharacter).*A_t);
    case 5 % symmetrized Itakura–Saito-type divergence. shen 2025 topology
        ratio = uniformCharacter./normGradientVec.^2;
        totalArea = sum(A_t);
        Average = sum(ratio.*A_t)/totalArea;
        basicEnergy = sum((ratio./Average + Average./ratio).*A_t)/totalArea;
    otherwise
        basicEnergy = sum ((vecnorm(gradientVec,2,2)-sqrt((k_c + curvature2Gradient)/8)).^2 .*A_t);
end
