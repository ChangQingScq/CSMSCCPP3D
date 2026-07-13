function [vCSM] = conformalSlitMappingGNK(vFlatten, freeboundsList, cornerIdx, longestBoundIdx, n, originOffset, alpha);
% Generalized Neumann Kernel method
%mesh3DWaveSurface
format long g
%
addpath bie fmm; 
tic
%%
% vertices of the out polygon: oriented counterclockwise
et = [];
etp = [];
for idx = 1:length(freeboundsList)
    [t_et, t_etp, ~] = splineReparamization(vFlatten(freeboundsList{idx}, :), n, cornerIdx{idx});
    if idx == longestBoundIdx
        et = [t_et(:,1) + sqrt(-1) * t_et(:,2);et];
        etp = [t_etp(:,1) + sqrt(-1) * t_etp(:,2);etp];
    else
        et = [et; t_et(:,1) + sqrt(-1) * t_et(:,2)];
        etp = [etp; t_etp(:,1) + sqrt(-1) * t_etp(:,2)];
    end
end

et  = et - originOffset;
% 
m = length(et)/n - 1; %the domain is of connectivity m+1
%
% The origin "0" must be in the domain.
%


%% 打印调试
% % Plot the origional domain
% colrr = ['k';'r';'b';'g';'m'];
% %colrr=jet(m+1);
% figure(1);
% hold on
% box off
% for k=1:m+1
%     crv = et((k-1)*n+1:k*n,1);
%     plot(real(crv),imag(crv),'Color',colrr(k),'LineWidth',4)
% end
% %plot(real(z1),imag(z1),'ko','MarkerSize',5)
% plot(0,0,'ko','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k')
% axis equal
% set( gca, 'XTick', [], 'YTick', [] );
% set (gca,'Visible','off');
% %axis([-75-real(originOffset) 75-real(originOffset) -50-imag(originOffset) 50-imag(originOffset)])
%%
%
%
%
% ===concentric circles slit mapping ===
% As in section 4.3 in the paper: 
% [N1] M.M.S. Nasser, Numerical conformal mapping of multiply connected
% regions onto the second, third and fourth  categories of Koebe’s
% canonical slit domains, J. Math. Anal. Appl. 382 (2011) 47–56.     
% Here, we choose: (theta_k=pi/2 for all k), since the slits are circular
% slits. You can change the values of thet_k to get different canonical
% domain
% 
thetk = ones(m + 1,1) * pi/2;
%
% The function A (see equation (4) in the above paper).
for k=1:m+1
    thet(1+(k-1)*n:k*n,1) = thetk(k);
end
A  = exp(1i*(pi/2-thet)).*et;
Ap = exp(1i*(pi/2-thet)).*etp;
%
% Compute the function gamma (the domain is bounded and the function gamma
% is given by (42a) in [N1].

if 1 == length(alpha)
    gam = -real(exp(1i*(pi/2-thet)).*clog(1-et/alpha));
else
    gam = -real(exp(1i*(pi/2-thet)).*clog(et));
end
%
% Solve the integral equation for "mu" and "h" (see equations (11) and (12)
% in [N1]. The integral equation is solved using the MATLAB function
% presented in:
% [N2] M.M.S. Nasser, Fast solution of boundary integral equations with the
% generalized Neumann kernel, Electronic Transactions on Numerical
% Analysis, Volume 44, pp. 189–229, 2015.  
%
[mun,h]=fbie(et,etp,A,gam,n,5,[],1e-13,1000);

%
% compute the constants h_k
for k=1:m+1
    hk(k,1)=mean(h(1+(k-1)*n:k*n));
end
c  =  exp(-hk(1));
%
% compute the boundary values of "f" (equation (10) in N1).
fet = (gam+h+1i*mun)./A;
% 
% compute the boundary values of the mapping function (equation (37) in N1)
%

if 1 == length(alpha)
    wet = c.*(1-et/alpha).*exp(et.*fet);
else
    wet = c.*et.*exp(et.*fet);
end
%
vFlattenC = vFlatten(:, 1) + sqrt(-1) * vFlatten(:, 2) - originOffset;
fv   = fcau(et,etp,fet, vFlattenC.');

if 1 == length(alpha)
    vCSM    = c.*(1-vFlattenC/alpha).'.*exp(vFlattenC.'.*fv);
else
    vCSM = c.*(vFlattenC).'.*exp((vFlattenC).'.*fv);
end

%% 一个不影响使用的bug，产生原因为柯西插值时n_bar内的点与n内的点过于接近引起的放大误差
IvalidNum = find(isnan(vCSM)==1)';
IvalidPt = vFlatten(IvalidNum,:);
NearestPtOnBoundaryNum = dsearchn([real(et),imag(et)],IvalidPt);
vCSM(IvalidNum) = wet(NearestPtOnBoundaryNum);
%%


