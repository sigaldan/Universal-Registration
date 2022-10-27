function [C_slave2master, angles_rad, cost_init, cost_post] = solve_wahba_svd(trk1, trk2, plott)
angles_rad = zeros(1,3);
if ~isempty(trk1) && ~isempty(trk2)
    num_sns_meas = length(trk1(:,1));
    
    r = trk2;
    b = trk1;
    
    r_iter = r;
    b_iter = b;
    if plott
        figure;
    end
    if plott
        figure(gcf); plot3(r_iter(:,1),r_iter(:,2),r_iter(:,3),'LineWidth',2);
        hold on;
        set(gca, 'ZDir','Reverse')
        plot3(b_iter(:,1),b_iter(:,2),b_iter(:,3),'g')
    end
    
    cost_init = 0;
    for i=1:num_sns_meas
        cost_init = cost_init + (norm(trk1(i,:)-trk2(i,:)))^2;
    end
    
    for j=1:1
        B = zeros(3);
        for i=1:num_sns_meas
            B = B + ((b_iter(i,:)') * (r_iter(i,:)));
        end
        %-----------------------------------------------------
        %         plott = 0;
        [U3 , ~, V3] = svd(B);
        C3 = U3 * diag([1 1 det(U3)*det(V3)]) * V3';
        C3 = U3 * diag([1 1 1]) * V3';
        
        r_iter_rot = (C3 * r_iter')';
        if plott
            plot3(r_iter_rot(:,1),r_iter_rot(:,2),r_iter_rot(:,3),'r')
        end
    end
    view(2)
    C_slave2master = C3;
    if plott
        close(gcf);
    end
    
else
    C_slave2master = eye(3);
end

cost_post = 0;
for i=1:num_sns_meas
    cost_post = cost_post + (norm(trk1(i,:)-(C_slave2master*trk2(i,:)')'))^2;
end

1;