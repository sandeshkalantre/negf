%NEGF calculation for I_J spectrum near the superconducting gap

%tolerance for calculation of surface Green's function
eps = 1e-4;

%0+ to get convergence to retarded Green's function
eta = 1e-8;

%energies in eV
t0 = 0.1;
Delta = 0.008;
kT = 0.002;
V = 0.0001;

phi = pi/2;

Delta1 = Delta;
alpha1 = [2*t0 Delta1; conj(Delta1) -2*t0];
beta1 = -t0* [1 0; 0 -1];

Delta2 = Delta * exp(1i*phi);
alpha2 = [2*t0+V Delta2; conj(Delta2) -2*t0+V];
beta2 = -t0* [1 0; 0 -1];

E = Delta*(-10:0.13:2);
I_QP = zeros(1,length(E));

for ii = 1:length(E)
    EE = E(ii);
    
    %calculation of g iteratively
    g1 = inv(alpha1);
    g1_last = inv(alpha1);
    err = 1;
    %calculate g1 for E = E
    for kk = 1:1000
        g1 = inv((EE + 1i*eta)*eye(2) - alpha1 - beta1'*g1*beta1);
        %for faster convergence
        err = sum(sum(abs(g1_last - g1)))/sum(sum(abs(g1) + abs(g1_last)));
        if  err < eps
            break
        end
        g1 = 0.5 * (g1 + g1_last);
        g1_last = g1;
    end
    
    g2 = inv(alpha2);
    g2_last = inv(alpha2);
    err = 1;
    %calculate g2 for E = E + V
    for kk = 1:1000
        g2 = inv((EE + V + 1i*eta)*eye(2) - alpha2 - beta2'*g2*beta2);
        %for faster convergence
        g2 = 0.5 * (g2 + g2_last);
        err = norm(g2 - g2_last,1)/norm(g2,1);
        if  err < eps
            break
        end
        g2_last = g2;
    end
    
    a1 = 1i * (g1 - g1');
    a2 = 1i * (g2 - g2');
    
    I_QP(ii) = 2 * a2(1,1) * a1(1,1) * (1.0/(1 + exp(EE/kT)) - 1.0/(1 + exp((EE+V)/kT)));
end

plot(E./Delta,I_QP,'LineWidth',2);
set(0,'DefaultTextInterpreter', 'latex');
set(gca,'Fontsize',[24]);
xlabel('$\frac{E}{\Delta}$','FontSize',24);
ylabel('$I_{QP}$','FontSize',24);
legend_handle = legend('$\phi = \frac{\pi}{2}, V = 0.001$');
set(legend_handle,'Interpreter','latex');

    
    
        