function [b1, b2] = pinBoundary(V, F)
%% A helper function to pin 2 points on the boundary
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%% Returns:
%%      b1: pinned boundary vertex
%%      b2: pinned boundary vertex

[B, ~] = findBoundary(V, F);
nB = length(B);

%% select the first and middle boundary points
b1 = B(1); b2 = B(round(nB/2));

end