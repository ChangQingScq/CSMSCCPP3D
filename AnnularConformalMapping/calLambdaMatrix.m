function [X, lambdaMatrix] = calLambdaMatrix(all_omega, outer_bd_id, inner_bd_id, v, f)

bdNum = length(all_omega)/2;

wedgeProductMatrix = zeros(bdNum * 2);
wedgeStarProductMatrix = zeros(bdNum * 2);
for idx = 1 : (bdNum * 2)
    for idx2 = 1 : (bdNum * 2)
        if idx == outer_bd_id || idx2 == outer_bd_id || ...
           idx == outer_bd_id + bdNum || ...
           idx2 == outer_bd_id + bdNum
            continue;
        end
        [energy1, energy2] = calWedgeProductEnergy(v, f, all_omega{idx}, all_omega{idx2});
        wedgeProductMatrix(idx, idx2) = energy1;
        wedgeStarProductMatrix(idx, idx2) = energy2;
    end
end

wedgeProductMatrix(:,outer_bd_id + bdNum) = [];
wedgeProductMatrix(outer_bd_id + bdNum,:) = [];
wedgeStarProductMatrix(:,outer_bd_id + bdNum) = [];
wedgeStarProductMatrix(outer_bd_id + bdNum,:) = [];

wedgeProductMatrix(:,outer_bd_id) = [];
wedgeProductMatrix(outer_bd_id,:) = [];
wedgeStarProductMatrix(:,outer_bd_id) = [];
wedgeStarProductMatrix(outer_bd_id,:) = [];

lambdaMatrix = wedgeProductMatrix\wedgeStarProductMatrix;


A = lambdaMatrix(length(all_omega)/2:end,1:length(all_omega)/2-1);
b = zeros(length(all_omega)/2, 1);
b(inner_bd_id) = 2*pi;
b(outer_bd_id) = [];
% b(outer_bd_id) = 2*pi;
% b(inner_bd_id) = [];
X = A \ b;

X = insertrows(X,0,outer_bd_id -1);

lambdaMatrix = insertrows(lambdaMatrix, zeros(1, size(lambdaMatrix, 2)), outer_bd_id-1);
lambdaMatrix = insertrows(lambdaMatrix, zeros(1, size(lambdaMatrix, 2)), outer_bd_id + length(all_omega)/2-1);
lambdaMatrix = insertrows(lambdaMatrix', zeros(size(lambdaMatrix, 1), 1)', outer_bd_id-1)';
lambdaMatrix = insertrows(lambdaMatrix', zeros(size(lambdaMatrix, 1), 1)', outer_bd_id + length(all_omega)/2-1)';



