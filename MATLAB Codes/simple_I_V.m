%NEGF calculation for a simple tunnel junnction between two normal metals
%tolerance for calculation of surface Green's function
eps = 1e-4;

%0+ to get convergence to retarded Green's function
eta = 1e-4;

%energies in eV
t0 = 0.1;
kT = 0.002;
mu = 2*t0;

V_vec = 0.001*(-1:0.1:1);
I_V = zeros(1,length(V_vec));

for jj = 1:length(V_vec)
    V = V_vec(jj);
    
    alpha1 = 2*t0 + V - mu;
    beta1 = -t0;
    
    alpha2 = 2*t0 - mu;
    beta2 = -t0;
    
    E_max = 0.001;
    E = E_max*(-1:0.1:1);
    I_E = zeros(1,length(E));
    a1_vec= zeros(1,length(E));
    a2_vec= zeros(1,length(E));
    
    %calculation of surface Green's function
    for ii = 1:length(E)
        EE = E(ii);

        %calculation of g iteratively
        g1 = inv(alpha1 + eps);
        g1_last = inv(alpha1 + eps);
        %calculate g1 for E = E
        for kk = 1:1000
            g1 = inv((EE + 1i*eta) - alpha1 - beta1'*g1*beta1);
            %for faster convergence
            err = sum(sum(abs(g1_last - g1)))/sum(sum(abs(g1) + abs(g1_last)));
            if  err < eps
                break
            end
            g1 = 0.5 * (g1 + g1_last);
            g1_last = g1;
        end

        g2 = inv(alpha2 + eps);
        g2_last = inv(alpha2 + eps);
        err = 1;
        %calculate g2 for E = E + V
        for kk = 1:1000
            g2 = inv((EE  + 1i*eta) - alpha2 - beta2'*g2*beta2);
            %for faster convergence
            err = sum(sum(abs(g2_last - g2)))/sum(sum(abs(g2) + abs(g2_last)));
            if  err < eps
                break
            end
            g2 = 0.5 * (g2 + g2_last);
            g2_last = g2;
        end

        a1 = 1i * (g1 - g1');
        a2 = 1i * (g2 - g2');
        
        a1_vec(ii) = a1;
        a2_vec(ii) = a2;
        
        I_E(ii) = 2 * a2(1,1) * a1(1,1) * (1.0/(1 + exp((EE-V)/kT)) - 1.0/(1 + exp((EE)/kT))); 
    end
    %plot(E,a1_vec,'LineWidth',2);
    %hold on;
    %plot(E,a2_vec,'LineWidth',2);
    %hold off;
    I_V(jj) = sum(I_E) * E_max * 0.1;
end
plot(V_vec,I_V,'LineWidth',2);

    
