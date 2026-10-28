function[energy1, energy2] = calWedgeProductEnergy(v, f, W1, W2)

energy1 = 0;
energy2 = 0;
for idx = 1:length(f)
    t_f = circshift([1 2 3], -1);
    t_f_shift = circshift([1 2 3], 1);
    [cot_A, cot_B, cot_C] = calPatchCot(v(f(idx,:), :));
    energy1 = energy1 + 0.5*det([W1(idx, t_f_shift) - W1(idx, t_f);W2(idx, t_f_shift) - W2(idx, t_f); 1, 1, 1]);
    energy2 = energy2 + 0.5*sum((W1(idx, t_f_shift) - W1(idx, t_f)).*(W2(idx, t_f_shift) - W2(idx, t_f)).* ...
              [cot_A, cot_B, cot_C]);
end


function[cot_A, cot_B, cot_C] = calPatchCot(v)
    A = v(1,:);
    B = v(2,:);
    C = v(3,:);
    
    AB = B - A;
    BC = C - B;
    CA = A - C;
    BA = A - B;
    CB = B - C;
    AC = C - A;
    % 计算边长
    a = norm(BC); % BC
    b = norm(CA); % CA
    c = norm(AB); % AB
    
    % 计算每个角的余弦值（使用点积公式）
    cos_A = -dot(AB, AC) / (c * b);
    cos_B = -dot(BC, BA) / (a * c);
    cos_C = -dot(CA, CB) / (b * a);

    cot_A = cotd(rad2deg(acos(cos_A)));
    cot_B = cotd(rad2deg(acos(cos_B)));
    cot_C = cotd(rad2deg(acos(cos_C)));

