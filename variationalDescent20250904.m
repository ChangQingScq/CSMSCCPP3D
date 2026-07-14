function[vertexEnergy, vertexEnergyGradient] = variationalDescent20250904(f, v, vCSM, lambda, FaceSFMMatrix, Tp, Bp, EsMode, k_c, vertex_triangles, ordered_ring, Avertex)
FV.faces = f;
FV.vertices = v;
[freeboundsList, ~] = getMeshBounds(f, v);

freeboundIdx = [];
for boundIdx = 1:length(freeboundsList)
    freeboundIdx = [freeboundIdx; freeboundsList{boundIdx}];
end

scalar = abs(vCSM);
%scalar = imag(vCSM); %debug
vertexOnBound = ismember([1:length(v)]',freeboundIdx);

%% 

vertexEnergy_Kw = zeros(length(v),1);
vertexEnergy_Kn = zeros(length(v),1);
vertexEnergy_Kg = zeros(length(v),1);
d_vertexEnergy_Kw = zeros(length(v),1);
d_vertexEnergy_Kn = zeros(length(v),1);
d_vertexEnergy_Kg = zeros(length(v),1);

for vIdx = 1:size(v, 1)
    f_temp = f(vertex_triangles{vIdx},:);
    [gradientVec, gradientTangent, normalizedNormVec, triangleCenter, A_t] = compute_gradient(v, f_temp, scalar');
    [d_g, d_gt, normalizedNormVec, triangleCenter, A_t] = compute_gradient_diff(v, f_temp, vIdx);
    norm_g = vecnorm(gradientVec, 2, 2);
    norm_gt = vecnorm(gradientTangent, 2, 2);
    unit_g = gradientVec ./ (norm_g + eps);
    %unit_gt = gradientTangent ./ (norm_gt + eps);
    d_norm_g = dot(d_g, unit_g, 2);
    %d_norm_gt = dot(d_gt, unit_gt, 2);
    d_unit_g = d_g ./ (norm_g + eps) - (gradientVec .* d_norm_g)./ (norm_g.^2 + eps);
    %d_unit_gt = d_gt ./ (norm_gt + eps) - (gradientTangent .* d_norm_gt)./ (norm_gt.^2 + eps);

    uComponent = dot(Tp(vertex_triangles{vIdx},:), unit_g, 2);
    vComponent = dot(Bp(vertex_triangles{vIdx},:), unit_g, 2);
    uDiffComponent = dot(Tp(vertex_triangles{vIdx},:), d_unit_g, 2);
    vDiffComponent = dot(Bp(vertex_triangles{vIdx},:), d_unit_g, 2);
    %% cal vertexEnergy_Kw & d_vertexEnergy_Kw
    if 0 ~= lambda(1)

        k_s = uComponent.^2.*FaceSFMMatrix(vertex_triangles{vIdx},1) + 2*uComponent.*vComponent.*FaceSFMMatrix(vertex_triangles{vIdx},2) + vComponent.^2.*FaceSFMMatrix(vertex_triangles{vIdx},3);
        k_s = -k_s;
        %dk_s = uDiffComponent.^2.*FaceSFMMatrix(vertex_triangles{vIdx},1) + 2*uDiffComponent.*vDiffComponent.*FaceSFMMatrix(vertex_triangles{vIdx},2) + vDiffComponent.^2.*FaceSFMMatrix(vertex_triangles{vIdx},3);    
        dk_s = 2*uDiffComponent.*uComponent.*FaceSFMMatrix(vertex_triangles{vIdx},1) + 2*uComponent.*vDiffComponent.*FaceSFMMatrix(vertex_triangles{vIdx},2) + 2*uDiffComponent.*vComponent.*FaceSFMMatrix(vertex_triangles{vIdx},2) + 2*vDiffComponent.*vComponent.*FaceSFMMatrix(vertex_triangles{vIdx},3);    
        dk_s = -dk_s;
        spacingScallopHeightMatchEnergy = 1; % Zou Qiang, computer aided design, 2014
        spacingUniformAntisymEnergy = 2; % Shen Changqing, 2025
        switch EsMode(1)
            case spacingScallopHeightMatchEnergy
                vertexEnergy_Kw(vIdx) = sum ( (norm_g-sqrt((k_s + k_c)/8)).^2 .*A_t );
    
                dKw = 2 * (norm_g - sqrt((k_s + k_c) / 8)) .* (d_norm_g - 0.5 * ((k_s + k_c) ./ 8).^0.5 .* dk_s);
                d_vertexEnergy_Kw(vIdx) = sum(dKw .* A_t);
            case spacingUniformAntisymEnergy
                vertexEnergy_Kw(vIdx) = sum ( ...
                    ((k_s + k_c)./8./norm_g.^2 + 8*norm_g.^2./(k_s + k_c)) .* A_t);
    
                dKw = 1/8*((dk_s ./ norm_g.^2) - 2*(k_s + k_c).*norm_g.^-3.*d_norm_g) + ...
                    8*((2*norm_g.*d_norm_g)./(k_s + k_c) - norm_g.^2.*(k_s + k_c).^-2.*dk_s);
                d_vertexEnergy_Kw(vIdx) = sum(dKw .* A_t);
            case 3
                vertexEnergy_Kw(vIdx) = sum ( ...
                    ((k_s + k_c)./8./norm_g.^2./0.3 + 8*norm_g.^2.*0.3./(k_s + k_c)) .* A_t);
    
                dKw = 1/8*((dk_s ./ norm_g.^2) - 2*(k_s + k_c).*norm_g.^-3.*d_norm_g)./0.3 + ...
                    8*((2*norm_g.*d_norm_g)./(k_s + k_c) - norm_g.^2.*(k_s + k_c).^-2.*dk_s).*0.3;
                d_vertexEnergy_Kw(vIdx) = sum(dKw .* A_t);
            otherwise
                vertexEnergy_Kw(vIdx) = sum ( (norm_g-sqrt((k_s + k_c)/8)).^2 .*A_t );
    
                dKw = 2 * (norm_g - sqrt((k_s + k_c) / 8)) .* (d_norm_g - 0.5 * ((k_s + k_c) ./ 8).^0.5 .* dk_s);
                d_vertexEnergy_Kw(vIdx) = sum(dKw .* A_t);
        end
    end
end

if 0 ~= lambda(2)
    [vertexEnergy_Kn, vertexEnergy_Kg, d_vertexEnergy_Kn, d_vertexEnergy_Kg] = getEk_dEK(f, v, vCSM, FaceSFMMatrix, Tp, Bp, vertex_triangles, ordered_ring, Avertex);
    vertexEnergy_Kg(freeboundIdx) = 0;  %边界上的Kg_diff不变
    %d_vertexEnergy_Kg(freeboundIdx) = 0;  %边界上的Kg_diff不变
end

%% cal vertexConformalEnergy
vertexConformalEnergy = zeros(length(v), 1);
ConformalEnergyGradient_uv = zeros(size(v, 1), 2);
if 0 ~= lambda(3)
    symmetricDirichletEnergy = 1; %对称迪利克雷能量
    symmetricConformalDeformEnergy = 2; %对称共形畸变能量

    for vIdx = 1:size(v, 1)
        f_pt = f(vertex_triangles{vIdx},:);
        for f_idx = 1:size(f_pt, 1) %点2是当前点
            f_pt(f_idx, :) = circshift(f_pt(f_idx,:), 2 - find(vIdx == f_pt(f_idx,:)));
        end
        E_sym = 0;
        dE_sym_du = 0;
        dE_sym_dv = 0;
        for fIdx = 1:size(f_pt, 1)
            f_temp = f_pt(fIdx, :);
            %oneFace = circshift(f_temp(facesIdx,:), 1 - find(vIdx == f_temp(facesIdx,:)));
            v_temp = v(f_temp, :);
            vCSM_temp = vCSM(f_temp);
            vCSM_temp = [real(vCSM_temp); imag(vCSM_temp)].';
    
            e1 = v_temp(2,:).' - v_temp(1,:).';
            e2 = v_temp(3,:).' - v_temp(1,:).';
            At = 0.5 * norm(cross(e1, e2));
    
            up = e1 / norm(e1);
            n = cross(e1, e2);
            n = n / norm(n);
            vp = cross(n, up);
    
            V2D = zeros(2, 3);
            V2D(:,1) = [0; 0];
            V2D(:,2) = [norm(e1); 0];
            vec = e2;
            V2D(:,3) = [dot(vec, up); dot(vec, vp)];
    
            E  = [V2D(:,2) - V2D(:,1), V2D(:,3) - V2D(:,1)];  % 2x2
            Ep = [vCSM_temp(2,:)' - vCSM_temp(1,:)', vCSM_temp(3,:)' - vCSM_temp(1,:)'];  % 2×2，每列是边向量
            J = Ep / E;  % 2×2
            G = J' * J;  % Gram 矩阵
    
            dEp_du = [1 0; 0 0];%vCSM_temp(:,2)是微扰点
            dJ_du = dEp_du / E;
            dG_du = dJ_du'*J + J'*dJ_du;
    
            dEp_dv = [0 0; 1 0];  % 2x2, 只扰动 Ep 的第一列的第二个分量
            dJ_dv = dEp_dv / E;
            dG_dv = dJ_dv' * J + J' * dJ_dv;
            % 检查奇异性
            if det(G) < 1e-50
                warning(sprintf('Degenerate triangle mapping. vIdx=%d, fIdx=%d',vIdx,fIdx));
                %E_sym = 10^60 - log(det(G)); %无穷大
                %continue;
            end
    
            switch EsMode(2)
                case symmetricDirichletEnergy
                    invG = inv(G);
                    JF2    = trace(G);
                    JinvF2 = trace(invG);
            
                    E_sym = E_sym + 0.5 * (JF2 + JinvF2) * At;
    
                    dE_sym_dG = 0.5 * (eye(2) - invG * invG);  % 2x2
                    dE_sym_du = dE_sym_du + trace(dE_sym_dG' * dG_du) * At;
                    dE_sym_dv = dE_sym_dv + trace(dE_sym_dG' * dG_dv) * At;
                case symmetricConformalDeformEnergy
                    % Step 3: Gram 矩阵 G = JᵀJ
                    [U, S, V] = svd(J);
                    sigma = diag(S);
                    sigma1 = sigma(1);
                    sigma2 = sigma(2);
            
                    E_sym = E_sym + (sigma1 / sigma2 + sigma2 / sigma1) * At;
    
                    dE_dSigma = diag([
                                1/sigma2 - sigma2 / sigma1^2;
                                1/sigma1 - sigma1 / sigma2^2
                            ]);
                    dE_dJ = U * dE_dSigma * V';
    
                    dE_sym_du = dE_sym_du + sum(sum(dE_dJ .* dJ_du)) * At;  % 等价于trace(dE_dJ' * dJ_du) * At
                    dE_sym_dv = dE_sym_dv + sum(sum(dE_dJ .* dJ_dv)) * At;  % 等价于trace(dE_dJ' * dJ_dv) * At
            end
        end
        ConformalEnergyGradient_uv(vIdx, :) = [dE_sym_du, dE_sym_dv]; % 梯度上升的反向
        vertexConformalEnergy(vIdx) = E_sym;
        %vertexConformalEnergy(vIdx) = log(E_sym); % log更容易收敛
    end
end
%%
vertexEnergy = vertexEnergy_Kw * lambda(1) + ...
               vertexEnergy_Kg * lambda(2) + ...
               vertexEnergy_Kn * lambda(2) + ...
               vertexConformalEnergy * lambda(3);
%%
vCSM_XY = [real(vCSM)', imag(vCSM)'];

if ~isempty(find(vecnorm(vCSM_XY, 2, 2) < eps))
    vCSM_XY(find(vecnorm(vCSM_XY, 2, 2) < eps), :) = 1e-9 * ConformalEnergyGradient_uv(find(vecnorm(vCSM_XY, 2, 2) < eps), :) ./ (eps + vecnorm(ConformalEnergyGradient_uv(find(vecnorm(vCSM_XY, 2, 2) < eps), :), 2, 2)); %保护过于靠近0的点vCSM_XY无方向，优先ConformalEnergyGradient_uv向量
end
if ~isempty(find(vecnorm(vCSM_XY, 2, 2) < eps))
    vCSM_XY(find(vecnorm(vCSM_XY, 2, 2) < eps), :) = [0, 1e-9]; %保护过于靠近0的点vCSM_XY无方向
end

norm_vCSM_XY = vCSM_XY./ (vecnorm(vCSM_XY,2,2) + eps);


maxNorm_dKw = max([abs(d_vertexEnergy_Kw) ; 1]);
maxNorm_dKg_dKn = max(max(abs([d_vertexEnergy_Kg, d_vertexEnergy_Kn])));
maxNorm_dKc = max(vecnorm(ConformalEnergyGradient_uv, 2, 2));

maxNorm_Kw = 1;
maxNorm_Kg_Kn = 1;
maxNorm_Kc = 1;

maxNorm_dKg_dKn = 1;%debug
maxNorm_dKw = 1;%debug

gradientEw_Ek = d_vertexEnergy_Kw ./ (maxNorm_dKw + eps) * lambda(1) + ...
                d_vertexEnergy_Kg ./ (maxNorm_dKg_dKn + eps) * lambda(2) + ...
                d_vertexEnergy_Kn ./ (maxNorm_dKg_dKn + eps) * lambda(2);
vertexEnergyGradient = [gradientEw_Ek, gradientEw_Ek] .* (norm_vCSM_XY) + ...
               ConformalEnergyGradient_uv ./ (maxNorm_dKc + eps) * lambda(3);

%vertexEnergyGradient = [0*gradientEw_Ek, gradientEw_Ek] .* (norm_vCSM_XY) + ...
%               ConformalEnergyGradient_uv ./ (maxNorm_dKc + eps) * lambda(3);%debug

vertexEnergyGradient = -vertexEnergyGradient;