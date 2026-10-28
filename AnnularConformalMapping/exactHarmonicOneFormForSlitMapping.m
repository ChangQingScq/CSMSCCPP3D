function f_0 = exactHarmonicOneFormForSlitMapping(v, f, outer_bd_id, inner_bd_id)

bd = meshboundaries(f);

A = cotLaplacian(v, f);

B = [];
g = [];
for idx = 1:length(bd)
    B = [B; bd{idx}];
    g = [g; bd{idx}*0+(inner_bd_id == idx)];
end

nV = size(v, 1);

I = setdiff(1:nV, B);


AII = A(I,I);
AIB = A(I,B);
ABB = A(B,B);

a = AII \ (I'*0-AIB*g);

f_0 = zeros(nV, 1);
f_0(I) = a;
f_0(B) = g;

