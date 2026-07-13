function [B, N] = findBoundary(v, f)
TR = triangulation(f, v(:,1), v(:,2), v(:,3));

B = freeBoundary(TR);
N = find(B(:,2) == B(1));
B = B(1:N,1);
