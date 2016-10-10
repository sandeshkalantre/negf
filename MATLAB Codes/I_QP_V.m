%NEGF calculation for I_J vs V
%tolerance for calculation of surface Green's function
clear all;
eps = 1e-4;

%0+ to get convergence to retarded Green's function
eta = 1e-6;

%energies in eV
t0 = 0.1;
Delta = 0.008;
kT = 0.0002;
phi = 0.3*pi;
mu = 0;

%since we want to see V/Delta 
V_vec = 5 * Delta *(-1:0.05:1);
I_V = zeros(1,length(V_vec));

for ii = 1:length(V_vec)
    V = V_vec(ii)

    Delta1 = Delta;
    alpha1 = [2*t0-mu Delta1; conj(Delta1) -2*t0+mu];
    beta1 = -t0* [1 0; 0 -1];

    Delta2 = Delta * exp(1i*phi);
    alpha2 = [2*t0-mu Delta2; conj(Delta2) -2*t0+mu];
    beta2 = -t0* [1 0; 0 -1];

    %assuming Fermi level at 2*t0
    E = 10*Delta*(-1:0.001:1);
    I_QP = zeros(1,length(E));
    a1_vec= zeros(1,length(E));
    a2_vec= zeros(1,length(E));

    for jj = 1:length(E)
        EE = E(jj);

        %calculation of g iteratively
        %g1 = inv(alpha1 + eps .* eye(2));
        %g1_last = inv(alpha1 + eps.* eye(2));
        if jj == 1
            g1 = inv(alpha1);
            g1_last = inv(alpha1);
        end
        %calculate g1 for E = E
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

        %g2 = inv(alpha2 + eps .* eye(2));
        %g2_last = inv(alpha2 + eps.* eye(2));
        if jj == 1
            g2 = inv(alpha2);
            g2_last = inv(alpha2);
        end
        err = 1;
        %calculate g2 for E = E + V
        for kk = 1:1000
            g2 = inv((EE + V + 1i*eta)*eye(2) - alpha2 - beta2'*g2*beta2);
            err = norm(g2 - g2_last,1)/norm(g2,1);
            if  err < eps
                break
            end
            %for faster convergence
            g2 = 0.5 * (g2 + g2_last);
            g2_last = g2;
        end

        a1 = 1i * (g1 - g1');
        a2 = 1i * (g2 - g2');
        
        a1_vec(jj) = trace(a1);
        a2_vec(jj) = trace(a2);
     
        I_QP(jj) = (a2(1,1) * a1(1,1) + a2(2,2) * a1(2,2)) * (1.0/(1 + exp((EE)/kT)) - 1.0/(1 + exp((EE+V)/kT)));
    end
    %plot(E,a1_vec,'LineWidth',2);
    %hold on;
    %plot(E,a2_vec,'LineWidth',2);
    %hold off;
    I_V(ii) = sum(I_QP);
end

plot(V_vec/Delta,I_V,'LineWidth',2);
set(0,'DefaultTextInterpreter', 'latex');
set(gca,'Fontsize',[12]);
xlabel('$\frac{V}{\Delta}$','FontSize',24);
ylabel('$I_{QP}$','FontSize',24);
legend_handle = legend({strcat('$\Delta = \ $ ',num2str(Delta),' $kT =\  $',num2str(kT))},'FontSize',24,'Location','Best');
set(legend_handle,'Interpreter','latex');
