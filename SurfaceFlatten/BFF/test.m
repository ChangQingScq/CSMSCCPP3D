clear
close all
clc

%%求解初始化
model = createpde;
shape1 = [1,0,0,20]';
shape2 = [1,-10,0,5]';
shape3 = [1,10,0,5]';

%几何边界定义详见 Constructive Solid Geometry (CSG)
gd = [shape1,shape2,shape3];
sf = 'shape1-shape2-shape3';
%sf = 'shape3';
ns = char('shape1','shape2','shape3');
ns = ns';
[dl, bt] = decsg(gd,sf,ns);%返回g n条边界上的一个点
pg = geometryFromEdges(model,dl);

%%网格划分
Mesh = generateMesh(model,'Hmax',2,'GeometricOrder','linear','Hedge',{[1:length(dl)],0.01});
%Mesh = generateMesh(model,'Hmax',1,'GeometricOrder','linear');

[p,e,t] = meshToPet(Mesh);%详见Mesh Data as [p,e,t] Triples

% pdemesh(Mesh,'ElementLabels','on','NodeLabels','on');%打印网格
pdemesh(Mesh,'NodeLabels','off');


p = p';
t = t';
z = sqrt(401- p(:,1).^2 - p(:,2).^2);
v = [p,real(z)];
f = t(:,1:3);


%[f, v] = read_obj('C:\Users\Shen_Changqing\Desktop\自由曲面3.obj');
pt = BFFAuto(v, f);
[B, ~] = findBoundary(v, f);
u = zeros(length(B), 1);
%u(1:90) =1.5;
u= u+log(3);
pt2 = BFFScale(v, f, u);% see https://peng00bo00.github.io/blog/2022/GAMES301-BFF/ ---Fixed a small bug

%myPlot(pt,'.');

%[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(v,f);

TR1 = triangulation(f, v(:,1), v(:,2), v(:,3));
TR2 = triangulation(f, pt(:,1), pt(:,2), pt(:,2)*0);
TR3 = triangulation(f, pt2(:,1), pt2(:,2), pt2(:,2)*0-5);

figure(1)
axis equal
trimesh(TR1)
figure(2)
hold on
axis equal
trimesh(TR2)
figure(2)
axis equal
trimesh(TR3)

tic;
map = annulus_conformal_map(p,f);
toc;
