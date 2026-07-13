function [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,face)

% compute_curvature - compute principal curvature directions and values
%
%   [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,face,options);
%
%   Umin is the direction of minimum curvature
%   Umax is the direction of maximum curvature
%   Cmin is the minimum curvature
%   Cmax is the maximum curvature
%   Cmean=(Cmin+Cmax)/2
%   Cgauss=Cmin*Cmax
%   Normal is the normal to the surface
%
%   options.curvature_smoothing controls the size of the ring used for
%       averaging the curvature tensor.
%
%   The algorithm is detailed in 
%       David Cohen-Steiner and Jean-Marie Morvan. 
%       Restricted Delaunay triangulations and normal cycle. 
%       In Proc. 19th Annual ACM Symposium on Computational Geometry, 
%       pages 237-246, 2003. 
%   and also in
%       Pierre Alliez, David Cohen-Steiner, Olivier Devillers, Bruno Le巚y, and Mathieu Desbrun. 
%       Anisotropic Polygonal Remeshing. 
%       ACM Transactions on Graphics, 2003. 
%       Note: SIGGRAPH '2003 Conference Proceedings
%
%   Copyright (c) 2007 Gabriel Peyre

orient = 1;


[vertex,face] = check_face_vertex(vertex,face);

n = size(vertex,2);
m = size(face,2);

% associate each edge to a pair of faces 
i = [face(1,:) face(2,:) face(3,:)]; 
j = [face(2,:) face(3,:) face(1,:)];
L = sub2ind([n, n], i, j); 
[a,b]=hist(L,unique(L)); 
Q = ismember(L,b(a>1)); 
[i1, j1] = ind2sub([n, n], L(Q)); 
W = [i' j']; 
R = ismember(W, [i1' j1'], 'rows'); 
J = find(R); 
[Q, K] = sortrows(W(R,:)); 
J = J(K); 
[~, ia, ~] = unique(Q, 'rows'); 
R = true(size(J)); 
R(ia) = false; 
J = J(R);
s = [1:m 1:m 1:m]; 
s(J) = 0; 
A = sparse(i,j,s,n,n); % if (i,j) are an edge at some face then A(i,j) is 
                    % this face index.
[i,j,~] = find(A); % direct link
X = [i j]; 
Y = [j i]; 
[~, I1, ~]=intersect(X, Y, 'rows'); 
I = X(I1, :); 
K1 = I(I(:,1) < I(:,2),:); % only directed edges 
K2 = [K1(:,2) K1(:,1)]; 
[~, ~, s1]= find(A(sub2ind(size(A),K1(:,1), K1(:,2)))); 
[~, ~, s2]= find(A(sub2ind(size(A),K2(:,1), K2(:,2))));
% links edge->faces 
E = [s1 s2]; 
i=K2(:,1); 
j=K2(:,2); 
ne = length(i); % number of directed edges
% normalized edge 
e = vertex(:,j) - vertex(:,i); 
d = sqrt(sum(e.^2,1)); 
d(d<eps) = 1; 
e = e ./ repmat(d,3,1); 
% avoid too large numerics 
d = d./mean(d);

% normals to faces
[~,normal] = compute_normal(vertex,face);

% inner product of normals
dp = sum( normal(:,E(:,1)) .* normal(:,E(:,2)), 1 );
% angle un-signed
beta = acos(clamp(dp,-1,1));
% sign
cp = crossp( normal(:,E(:,1))', normal(:,E(:,2))' )';
si = orient * sign( sum( cp.*e,1 ) );
% angle signed
beta = beta .* si;
% tensors
T = zeros(3,3,ne);
for x=1:3
    for y=1:x
        T(x,y,:) = reshape( e(x,:).*e(y,:), 1,1,ne );
        T(y,x,:) = T(x,y,:);
    end
end
T = T.*repmat( reshape(d.*beta,1,1,ne), [3,3,1] );

% do pooling on vertices
Tv = zeros(3,3,n);
w = zeros(1,1,n);
for k=1:ne
%    progressbar(k,ne);
    Tv(:,:,i(k)) = Tv(:,:,i(k)) + T(:,:,k);
    Tv(:,:,j(k)) = Tv(:,:,j(k)) + T(:,:,k);
    w(:,:,i(k)) = w(:,:,i(k)) + 1;
    w(:,:,j(k)) = w(:,:,j(k)) + 1;
end
w(w<eps) = 1;
Tv = Tv./repmat(w,[3,3,1]);


% extract eigenvectors and eigenvalues
U = zeros(3,3,n);
D = zeros(3,n);
for k=1:n
    [u,d] = eig(Tv(:,:,k));
    d = real(diag(d));
    % sort acording to [norma,min curv, max curv]
    [~,I] = sort(abs(d));    
    D(:,k) = d(I);
    U(:,:,k) = real(u(:,I));
end

Umin = squeeze(U(:,3,:));
Umax = squeeze(U(:,2,:));
Cmin = D(2,:)';
Cmax = D(3,:)';
Normal = squeeze(U(:,1,:));
Cmean = (Cmin+Cmax)/2;
Cgauss = Cmin.*Cmax;

% enforce than min<max
I = find(Cmin>Cmax);
Cmin1 = Cmin; Umin1 = Umin;
Cmin(I) = Cmax(I); Cmax(I) = Cmin1(I);
Umin(:,I) = Umax(:,I); Umax(:,I) = Umin1(:,I);

% try to re-orient the normals
normal = compute_normal(vertex,face);
s = sign( sum(Normal.*normal,1) ); 
Normal = Normal .* repmat(s, 3,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = crossp(x,y)
% x and y are (m,3) dimensional
z = x;
z(:,1) = x(:,2).*y(:,3) - x(:,3).*y(:,2);
z(:,2) = x(:,3).*y(:,1) - x(:,1).*y(:,3);
z(:,3) = x(:,1).*y(:,2) - x(:,2).*y(:,1);
end

function y = clamp(x,a,b)

% clamp - clamp a value
%
%   y = clamp(x,a,b);
%
% Default is [a,b]=[0,1].
%
%   Copyright (c) 2004 Gabriel Peyr

if nargin<2
    a = 0;
end
if nargin<3
    b = 1;
end

y = max(x,a);
y = min(y,b);
end
