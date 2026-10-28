function[lambda] = calLambda(v, f, w1, w2)

energy1 = 0;
energy2 = 0;
for idx = 1:length(f)
    t_f = circshift(f(idx, :), -1);
    t_f_shift = circshift(f(idx, :), 1);
    [cot_A, cot_B, cot_C] = calPatchCot(v(f(idx,:), :));
    energy1 = energy1 + 0.5*det([w1(t_f_shift)' - w1(t_f)';w2(t_f_shift)' - w2(t_f)'; 1, 1, 1]);
%     energy2 = energy2 + 0.5*sum((w1(t_f_shift)' - w1(t_f)').*(w2(t_f_shift)' - w2(t_f)').* ...
%                         [W(t_f(2), t_f(3)),W(t_f(3), t_f(1)),W(t_f(1), t_f(2))]);
    energy2 = energy2 + 0.5*sum((w1(t_f_shift)' - w1(t_f)').*(w2(t_f_shift)' - w2(t_f)').* ...
                        [cot_A, cot_B, cot_C]);
end

lambda= energy1 ./ energy2;

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
    cos_A = dot(AB, AC) / (c * b);
    cos_B = dot(BC, BA) / (a * c);
    cos_C = dot(CA, CB) / (b * a);

    cot_A = cotd(rad2deg(acos(cos_A)));
    cot_B = cotd(rad2deg(acos(cos_B)));
    cot_C = cotd(rad2deg(acos(cos_C)));

