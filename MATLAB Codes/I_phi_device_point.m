%calculation of the Josephson current vs phase directly using the current
%operator
%device point added

%tolerance for calculation of surface Green's function
eps = 1e-4;

%0+ to get convergence to retarded Green's function
eta = 1e-8;

%energies in eV
t0 = 0.1;
Delta = 0.008;
kT = 0.002;

phi_vec = 2*pi*(0:0.05:1);
I_phi = zeros(2,length(phi_vec));

%tunelling
M0 = 0.1;
M = M0 .* [1 0; 0 1];

%device Hamiltonian
alpha = [2*t0 0; 0 -2*t0];
beta1 = -t0* [1 0; 0 -1];
H_D = alpha; 


for ii = 1:length(phi_vec)
    phi = phi_vec(ii);

    Delta1 = Delta;
    alpha1 = [2*t0 Delta1; conj(Delta1) -2*t0];
    beta1 = -t0* [1 0; 0 -1];

    Delta2 = Delta * exp(1i*phi);
    alpha2 = [2*t0 Delta2; conj(Delta2) -2*t0];
    beta2 = -t0* [1 0; 0 -1];

    E = 4 * t0 *(0:0.01:1);
    I_E = zeros(2,length(E));
    for jj = 1:length(E)
        EE = E(jj);
    
        %calculation of g iteratively
        if jj == 1
            g1 = inv(alpha1);
            g1_last = inv(alpha1);
        end
        err = 1;
        for kk = 1:1000
            g1 = inv((EE + 1i*eta)*eye(2) - alpha1 - beta1'*g1*beta1);
            err = norm(g1 - g1_last,1)/norm(g1,1);
            if  err < eps
                break
            end
            %for faster convergence
            g1 = 0.5 * (g1 + g1_last);
            g1_last = g1;
        end
    
        if jj == 1
            g2 = inv(alpha2);
            g2_last = inv(alpha2);
        end
        err = 1;
        for kk = 1:1000
            g2 = inv((EE + 1i*eta)*eye(2) - alpha2 - beta2'*g2*beta2);
            err = norm(g2 - g2_last,1)/norm(g2,1);
            if  err < eps
                break
            end
            %for faster convergence
            g2 = 0.5 * (g2 + g2_last);
            g2_last = g2;
        end
        
        g_D = inv((EE + 1i*eta) .* eye(2) - H_D);
        G_D = ...
        inv([inv(g1) -M zeros(2,2);-M' inv(g_D) -M;zeros(2,2) -M' inv(g2)]);
        sigma1_corr = (1.0/(1 + exp(EE/kT))) * beta1' * (g1' - g1) * beta1;
        sigma2_corr = (1.0/(1 + exp(EE/kT))) * beta2 * (g2' - g2) * beta2';
        Sigma_D_corr = ...
        [sigma1_corr zeros(2,2) zeros(2,2); 
         zeros(2,2) zeros(2,2) zeros(2,2); 
         zeros(2,2) zeros(2,2) sigma2_corr];
     
        G_corr = G_D * Sigma_D_corr * G_D';
        
        I_op = M * G_corr(3:4,1:2) - G_corr(1:2,3:4) * M;
        I_E(1,jj) = I_op(1,1) - I_op(2,2);    
        
        %calculation without device point, symbols resued, be careful in
        %further use
        G_D = inv([inv(g1) -M;-M' inv(g2)]);
        sigma1_corr = (1.0/(1 + exp(EE/kT))) * beta1' * (g1' - g1) * beta1;
        sigma2_corr = (1.0/(1 + exp(EE/kT))) * beta2 * (g2' - g2) * beta2';
        Sigma_D_corr = [sigma1_corr zeros(2,2);zeros(2,2) sigma2_corr];
        G_corr = G_D * Sigma_D_corr * G_D';
        
        I_op = M * G_corr(3:4,1:2) - G_corr(1:2,3:4) * M;
        I_E(2,jj) = I_op(1,1) - I_op(2,2);
        
    end
    I_phi(1,ii) = sum(I_E(1,:));
    I_phi(2,ii) = sum(I_E(2,:));
end

plot(phi_vec/(2*pi),I_phi(1,:),'LineWidth',2)
hold on;
%plot(phi_vec/(2*pi),I_phi(2,:),'LineWidth',2)
hold off;
set(0,'DefaultTextInterpreter', 'latex');
set(gca,'Fontsize',[16]);
xlabel('$\frac{\phi}{2\pi}$','FontSize',16);
ylabel('$I_J$ (arb units)','FontSize',16);
title(strcat('$M = \ $ ',num2str(M0)));
legend_handle = legend({'With Device Point','Without Device Point'},'Location','NorthEast');
set(legend_handle,'Interpreter','latex');
        