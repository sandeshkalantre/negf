%eigenvalues for tight binding with potential model
t0 = 0.1;
V = 0.0;

%number of lattice points
N = 1000;

%Hamiltonian
H = zeros(N,N);
for ii = 1:N
    H(ii,ii) = 2*t0 + V;
    if ii + 1 < N + 1
        H(ii,ii + 1) = -t0;
    end
    if ii - 1 > 0
        H(ii,ii - 1) = -t0;
    end
end

E_vec = sort(eig(H));
figure
subplot(2,1,1)
plot(E_vec,'LineWidth',2);
legend({strcat('t0 = ',num2str(t0),' V = ',num2str(V))},'FontSize',16,'Location','Best');
t0 = 0.1;
V = 1;

%number of lattice points
N = 1000;

%Hamiltonian
H = zeros(N,N);
for ii = 1:N
    H(ii,ii) = 2*t0 + V;
    if ii + 1 < N + 1
        H(ii,ii + 1) = -t0;
    end
    if ii - 1 > 0
        H(ii,ii - 1) = -t0;
    end
end

E_vec = sort(eig(H));

subplot(2,1,2)
plot(E_vec,'LineWidth',2);
legend({strcat('t0 = ',num2str(t0),' V = ',num2str(V))},'FontSize',16,'Location','Best');




