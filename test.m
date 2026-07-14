clear
clc
close all

%%
%matlab mesh-generation
%see https://www.mathworks.com/help/pde/ug/create-geometry-at-the-command-line.html?searchHighlight=Constructive%20Solid%20Geometry%20%28CSG%29&s_tid=srchtitle_support_results_1_Constructive%20Solid%20Geometry%20%2528CSG%2529

model = createpde;
shape1 = [1,0,0,10,10,0]';
shape2 = [1,7,0,1,0,0]';
shape3 = [1,-7,0,1,0,0]';
shape4 = [1, 0, 0, 1, 0, 0]';
gd = [shape1,shape2,shape3,shape4];
sf = 'shape1-shape2-shape3-shape4';
%sf = 'shape3';
ns = char('shape1','shape2','shape3','shape4');
ns = ns';
[dl, bt] =  decsg(gd,sf,ns);
pg = geometryFromEdges(model,dl);

Mesh = generateMesh(model,'Hmax', 1,'GeometricOrder','linear',"Hedge",{[1:16], 0.3});%1 0.3 v1048 f2475

[p,e,t] = meshToPet(Mesh); % see Mesh Data as [p,e,t] Triples
p = p';

%z = cos(sqrt(p(:,1).^2 + p(:,2).^2)/20) *0;
z = sin(pi/20*p(:,1)).*p(:,2)/10*5; %wave height
f = t(1:3,:)';
v = [p,z];




%% Flatten
flattenMode = "BaryCenter_BFF";

if "LSCM" == flattenMode
    vFlatten = LSCM(v, f);
    vFlatten(:,2) = -vFlatten(:,2);
elseif "BFF" == flattenMode
    [fout, vout, freeboundsList, longestBoundIdx] = holeFiller(f, v, "LSM");%maybe very slow
    vFlatten = BFFAuto(vout, fout);
    vFlatten = vFlatten(1:length(v), :);
    vFlatten = [vFlatten(:,2), vFlatten(:,1)];
elseif "BaryCenter_BFF" == flattenMode
    [fout, vout, freeboundsList, longestBoundIdx] = holeFiller(f, v, "BaryCenter");
    vFlatten = BFFAuto(vout, fout);
    vFlatten = vFlatten(1:length(v), :);
    vFlatten = [vFlatten(:,2), vFlatten(:,1)];
elseif "None" == flattenMode
    vFlatten = v(:,1:2);
else
    fprintf("invalid flatten mode:%s\n",flattenMode);
end

cornerIdx = cornerIdentify(v, freeboundsList);

[freeboundsList, longestBoundIdx] = getMeshBounds(f, v);
for idx = 1:length(freeboundsList)
        freeboundsList{idx} = flip(freeboundsList{idx});
end

%% Conformal slit mapping

ModeCSM = "GNKDISK";
vCSM = [];

if ModeCSM == "GNKDISK"
    n = 2^13; % a common multiple of cornerIdx
    originOffset = 3+1i*3; %originOffset should on S
    alpha = [];
    vCSM = conformalSlitMappingGNK(vFlatten, freeboundsList, cornerIdx, longestBoundIdx, n, originOffset, alpha);
elseif ModeCSM == "GNKANNULAR"
    n = 2^13; % a common multiple of cornerIdx
    originOffset = 3+3i;%originOffset should on S
    alpha = 0;%alpha should in a hole
    alpha = alpha-originOffset;
    vCSM = conformalSlitMappingGNK(vFlatten, freeboundsList, cornerIdx, longestBoundIdx, n, originOffset, alpha);
elseif ModeCSM == "H1FORM"
    inner_bd_id = 4;
    [tau] = annularSlitMapping(v, f, inner_bd_id);
    vCSM = exp(tau).';
end

vCSM = [real(vCSM)', imag(vCSM)'];

for idx = 1:length(freeboundsList)
    bdPt = freeboundsList{idx};
    radOnBd = vecnorm(vCSM(bdPt, :), 2, 2);
    avgRadOnBd = mean(radOnBd);
    vCSM(bdPt, :) = vCSM(bdPt, :) ./ (radOnBd + eps) * avgRadOnBd;
end

vCSM = (vCSM(:, 1) + 1i * vCSM(:, 2)).';


%%
TR1 = triangulation(f, real(vCSM)', imag(vCSM)'); %conformal slit mapping mesh
TR2 = triangulation(f, v(:,1), v(:,2), v(:,3)); %3D surface mesh
TR3 = triangulation(f, vFlatten(:,1), vFlatten(:,2)); %flattened surface mesh

figure(1)
hold on
axis equal
triplot(TR1,'color',[0.9 0.9 0.9])

figure(2)
hold on
axis equal
trimesh(TR2)

figure(3)
hold on
axis equal
triplot(TR3)

ptOnBound = [];
for idx = 1:length(freeboundsList)
    ptOnBound = [ptOnBound; freeboundsList{idx}];
end
ptNotOnBound = setdiff((1:length(vCSM))',ptOnBound);

%% Functional energy calculation

scalar = abs(vCSM); %test scalar

FV.faces = f;
FV.vertices = v;
[FaceSFM, Bp, Tp] = GetSecondBasicForm( FV );
FaceSFMMatrix = zeros(length(FaceSFM),3);
for idx = 1:length(FaceSFM)
    curFaceSFM = cell2mat(FaceSFM(idx));
    FaceSFMMatrix(idx,:) = [curFaceSFM(1,1), curFaceSFM(1,2), curFaceSFM(2,2)];
end

[~, ~, normalizedNormVec, ~, ~] = compute_gradient(v, f, scalar');
[~, Avertex, ~, ~, ~] = CalcVertexNormals(FV, normalizedNormVec);
allE = [];
lambda = [1 0 0]; % [Ew Ek Ec] %Ec conformal energy
EsMode = [2, 2];
k_c = 1;
vertexEnergyGradient = zeros(length(vCSM),2);
adaptionGrad = vertexEnergyGradient;

scalar = abs(vCSM);
[curvature2Gradient, gradientVec, triangleCenter, A_t] = getCurvature2Gradient(v, f, scalar', Bp, Tp, FaceSFMMatrix);
energyW = calBasicEnergy(gradientVec, k_c, curvature2Gradient, A_t, EsMode(1));
energyK = getEk(v, f, scalar, Bp, Tp, FaceSFMMatrix);

alpha = 0.1;
E = energyW + alpha*energyK
