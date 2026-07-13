function [EFG,u,v] = compute_second_basic_form(vertex, face)
%%see DOI: 10.1109/TDPVT.2004.1335277

    [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex, face);

    Normal = Normal';

    node0 = vertex(face(:,1), :);
    node1 = vertex(face(:,2), :);
    node2 = vertex(face(:,3), :);
    e0 = node2 - node1;
    e1 = node0 - node2;
    e2 = node1 - node0;

    n0 = Normal(face(:,1), :);
    n1 = Normal(face(:,2), :);
    n2 = Normal(face(:,3), :);

    normalizedNormVec = cross(node1- node0, node2- node0)./vecnorm(cross(node1- node0, node2- node0),2,2); 
    u = e2./vecnorm(e2,2,2);
    v = cross(normalizedNormVec, u)./vecnorm(cross(normalizedNormVec, u),2,2);

    EFG = [];
    for idx = 1:size(u,1)
        A = [dot(e0(idx,:),u(idx,:)), dot(e0(idx,:),v(idx,:)), 0
             0, dot(e0(idx,:),u(idx,:)), dot(e0(idx,:),v(idx,:))
             dot(e1(idx,:),u(idx,:)), dot(e1(idx,:),v(idx,:)), 0
             0, dot(e1(idx,:),u(idx,:)), dot(e1(idx,:),v(idx,:))
             dot(e2(idx,:),u(idx,:)), dot(e2(idx,:),v(idx,:)), 0
             0, dot(e2(idx,:),u(idx,:)), dot(e2(idx,:),v(idx,:))
             ];

        %%

        B = [dot((n2(idx,:)-n1(idx,:)),u(idx,:)); 
             dot((n2(idx,:)-n1(idx,:)),v(idx,:))
             dot((n0(idx,:)-n2(idx,:)),u(idx,:))
             dot((n0(idx,:)-n2(idx,:)),v(idx,:))
             dot((n1(idx,:)-n0(idx,:)),u(idx,:))
             dot((n1(idx,:)-n0(idx,:)),v(idx,:))];
        
        EFG = [EFG;linsolve(A,B)'];
    end

end
