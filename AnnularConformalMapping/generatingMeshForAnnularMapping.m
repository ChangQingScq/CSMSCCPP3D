function[v, f] = generatingMeshForAnnularMapping()
model = createpde;
shape1 = [1,0,0,20]';
shape2 = [1,-10,0,5]';shape2 = [1,0,0,5]';
%shape3 = [1,10,0,5]';

%几何边界定义详见 Constructive Solid Geometry (CSG)
gd = [shape1,shape2];
sf = 'shape1-shape2';
%sf = 'shape3';
ns = char('shape1','shape2');
ns = ns';
[dl, bt] = decsg(gd,sf,ns);
pg = geometryFromEdges(model,dl);

%%网格划分
Mesh = generateMesh(model,'Hmax',1,'GeometricOrder','linear','Hedge',{[1:size(dl,2)],1});
%Mesh = generateMesh(model,'Hmax',1,'GeometricOrder','linear');

[p,e,t] = meshToPet(Mesh);%详见Mesh Data as [p,e,t] Triples

% pdemesh(Mesh,'ElementLabels','on','NodeLabels','on');%打印网格
pdemesh(Mesh,'NodeLabels','off');

p = p';
t = t';
z = sqrt(401- p(:,1).^2 - p(:,2).^2);%z = sqrt(401- p(:,1).^2 - p(:,2).^2);
v = [p,real(z)];
f = t(:,1:3);