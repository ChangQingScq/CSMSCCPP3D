function[tau] = calTau(all_omega, outer_bd_id, inner_bd_id, v, f);

[X, lambdaMatrix] = calLambdaMatrix(all_omega, outer_bd_id, inner_bd_id, v, f);

omega = all_omega{1} * 0;
for idx_X = 1:length(X)
    omega = omega + all_omega{idx_X} * X(idx_X);
    for idx_lambda = 1:length(lambdaMatrix)
        omega = omega + 1i * X(idx_X) * lambdaMatrix(idx_lambda, idx_X) * all_omega{idx_lambda};
    end
end

%% calculate d_tau
d_tau = zeros(length(f)*3,3);
for idx_f = 1:length(f)
    d_tau((idx_f-1)*3+1:(idx_f-1)*3+3,:) = [f(idx_f,:)
                                            circshift(f(idx_f,:),-1)
                                            circshift(omega(idx_f,:),-1)-omega(idx_f,:)]';
end

swapIdx = find(d_tau(:,1)>d_tau(:,2));
d_tau(swapIdx,:) = [d_tau(swapIdx,2),d_tau(swapIdx,1),-d_tau(swapIdx,3)]; 

[~,ia,~] = unique(d_tau(:,1:2),'rows');
d_tau = d_tau(ia,:);

%% path integret
tau = zeros(length(v),1);
tau_on_v_cal = tau + 1;
startPt = d_tau(1);
tau_on_v_cal(startPt) = 0;
while sum(tau_on_v_cal) ~= 0
    t_startPt = [];
    for idx = 1:length(startPt)
        connect_edge = find(startPt(idx) == d_tau(:,1));
        cur_d_imag_tau = d_tau(connect_edge, :);
        [connect_pt, ia, ~] = intersect(cur_d_imag_tau(:, 2), find(tau_on_v_cal == 1));
        tau(connect_pt) = (tau(startPt(idx)) + cur_d_imag_tau(ia,3));
        tau_on_v_cal(connect_pt) = 0;
        t_startPt = [t_startPt;connect_pt];

        connect_edge = find(startPt(idx) == d_tau(:,2));
        cur_d_imag_tau = d_tau(connect_edge, :);
        [connect_pt, ia, ~] = intersect(cur_d_imag_tau(:, 1), find(tau_on_v_cal == 1));

        tau(connect_pt) = (tau(startPt(idx)) - cur_d_imag_tau(ia,3));
        tau_on_v_cal(connect_pt) = 0;
        t_startPt = [t_startPt;connect_pt];
    end

    startPt = t_startPt;
end


