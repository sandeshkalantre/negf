%NEGF calculation for I_J vs V
%tolerance for calculation of surface Green's function
clear all;
eps = 1e-4;

%0+ to get convergence to retarded Green's function
eta = 1e-6;

%energies in eV
t0 = 0.1;
Delta = 0.003;
fac=1;
kT = 0.0001;
phi = 0.2*pi;
V_vec1 = Delta*(0:0.1:4);
I_V1 = zeros(1,length(V_vec1));
mu =0;

for jj = 1:length(V_vec1)
    V = V_vec1(jj);

    Delta1 = Delta;
    alpha1 = [2*t0-mu Delta1; conj(Delta1) -2*t0+mu];
    beta1 = t0* [1 0; 0 -1];

    Delta2 = fac * Delta * exp(1i*phi);
    alpha2 = [2*t0-mu Delta2; conj(Delta2) -2*t0+mu];
    beta2 = t0* [1 0; 0 -1];

    E = Delta*(-10:0.01:10);
    I_QP = zeros(1,length(E));

    for ii = 1:length(E)
        EE = E(ii);
        %calculation of g iteratively
        if ii==1
        g1 = inv(alpha1);
        g1_last = inv(alpha1);
        end
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
        if ii==1
        g2 = inv(alpha2);
        g2_last = inv(alpha2);
        end
        err = 1;
        %calculate g2 for E = E + V
        for kk = 1:1000
            g2 = inv((EE +V + 1i*eta)*eye(2) - alpha2 - beta2*g2*beta2');
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
        if jj==30
        DOS1(ii)=trace(a1);
        DOS2(ii)=trace(a2);
        end
 I_QP(ii) = (a2(1,1) * a1(1,1)+a2(2,2)*a1(2,2)) * (1.0/(1 + exp((EE)/kT)) - 1.0/(1 + exp((EE+V)/kT)));
    end
     I_V1(jj) = sum(I_QP);
end
V_vec2 = Delta*(0:-0.1:-4);
I_V2 = zeros(1,length(V_vec2));
mu =0;
for jj = 1:length(V_vec2)
    V = V_vec2(jj);

    Delta1 = Delta;
    alpha1 = [2*t0-mu Delta1; conj(Delta1) -2*t0+mu];
    beta1 = t0* [1 0; 0 -1];

    Delta2 = fac * Delta * exp(1i*phi);
    alpha2 = [2*t0-mu Delta2; conj(Delta2) -2*t0+mu];
    beta2 = t0* [1 0; 0 -1];

    E = Delta*(-10:0.01:10);
    I_QP = zeros(1,length(E));

    for ii = 1:length(E)
        EE = E(ii);
        %calculation of g iteratively
        if ii==1
        g1 = inv(alpha1);
        g1_last = inv(alpha1);
        end
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
        if ii==1
        g2 = inv(alpha2);
        g2_last = inv(alpha2);
        end
        err = 1;
        %calculate g2 for E = E + V
        for kk = 1:1000
            g2 = inv((EE +V + 1i*eta)*eye(2) - alpha2 - beta2*g2*beta2');
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
        if jj==30
        DOS1(ii)=trace(a1);
        DOS2(ii)=trace(a2);
        end
 I_QP(ii) = (a2(1,1) * a1(1,1)+a2(2,2)*a1(2,2)) * (1.0/(1 + exp((EE)/kT)) - 1.0/(1 + exp((EE+V)/kT)));
    end
     I_V2(jj) = sum(I_QP);
end
plot(E/Delta,DOS1,E/Delta,DOS2,'LineWidth',2);
xlim([-10 10]);
set(0,'DefaultTextInterpreter', 'latex');
set(gca,'Fontsize',[24]);
xlabel('$\frac{E}{\Delta}$','FontSize',24);
ylabel('$DOS$','FontSize',24);
legend_handle = legend('$\phi = \frac{\pi}{2}$');
set(legend_handle,'Interpreter','latex');
figure (2)
plot(V_vec1/Delta,I_V1,V_vec2/Delta,I_V2,'LineWidth',2);
set(0,'DefaultTextInterpreter', 'latex');
set(gca,'Fontsize',[24]);
xlabel('$\frac{V}{\Delta}$','FontSize',24);
ylabel('$I_{QP}$','FontSize',24);
legend_handle = legend('$\phi = \frac{\pi}{2}$');
set(legend_handle,'Interpreter','latex');
