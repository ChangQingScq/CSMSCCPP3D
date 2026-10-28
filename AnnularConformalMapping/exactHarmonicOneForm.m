function f_0 = exactHarmonicOneForm(v, f)

bd = meshboundaries(f); % get all boundaries
nb = length(bd); % number of boundaries

% determine which boundary is the outer one
bd1 = bd{1};
bd2 = bd{2};
if size(v,2) == 2
    % use the Pick's formula
    area1 = abs(sum(v(bd1,1).*v(bd1([2:length(bd1),1]),2)) - ...
        sum(v(bd1([2:length(bd1),1]),1).*v(bd1,2)));
    area2 = abs(sum(v(bd2,1).*v(bd2([2:length(bd2),1]),2)) - ...
        sum(v(bd2([2:length(bd2),1]),1).*v(bd2,2)));
    if area1 > area2
        outer_bd_id = 1;
        inner_bd_id = 2;
    else
        outer_bd_id = 2;
        inner_bd_id = 1;
    end        
else
    % pick the longer boundary as the outer boundary
    bd1_length = sum(sqrt(sum((v(bd1,:) - v(bd1([2:end,1]),:)).^2,2)));
    bd2_length = sum(sqrt(sum((v(bd2,:) - v(bd2([2:end,1]),:)).^2,2)));
    if bd1_length > bd2_length
        outer_bd_id = 1;
        inner_bd_id = 2;
    else 
        outer_bd_id = 2;
        inner_bd_id = 1;
    end
end

A = cotLaplacian(v, f);
B = [bd{outer_bd_id};bd{inner_bd_id}];
g = [bd{outer_bd_id}*0+1;bd{inner_bd_id}*0];
nV = size(v, 1);

I = setdiff(1:nV, B);


AII = A(I,I);
AIB = A(I,B);
ABB = A(B,B);

a = AII \ (I'*0-AIB*g);

f_0 = zeros(nV, 1);
f_0(I) = a;
f_0(B) = g;

