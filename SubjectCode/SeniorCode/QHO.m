%% Contents
hbar=6.58E-16; % [eV*s] 约化普朗克常数
c_cm=2.297E10; % [cm/s] 光速
c_nm=2.297E17; % [nm/s] 光速
amu2eV=931.5E6; % [eV] convert amu to eV [E=m*c^2]
%% 参数
N=9;%基态的数量
m_amu=6; % [amu] oscillator mass
m=m_amu*amu2eV/(c_nm^2); % [eV*s^2/nm^2] oscillator mass
f_inv_cm=300; % [cm^{-1}] oscillator frequency (f=c/lamda)
f=f_inv_cm*c_cm; % [Hz] oscillator frequency
w=2*pi*f; % [rad/s] oscillator angular frequency

nt = 75; % [] number of time points
TimeSpan = 3; % [T0] calculation time span in units of T0
%% 计算
% 构造QHO的算符
a=diag(sqrt(1:N-1),1); % 下降算符
ad=a.'; % 上升算符
H=hbar*w*(ad*a+0.5*eye(N)); % [eV] Hamiltonian (energy basis)
X=sqrt(hbar/(2*m*w))*(ad+a); % position operator (energy basis)
% UEQ - eigenvectors of the X operator (energy basis)
[UEX,Xeigval]=eig(X);
UXE=UEX'; % transformation matrix from energy to state basis
% To transform psi_E in energy basis to psi_X in position basis,use
% psi_E=UQE*psi_X

E0 = hbar*w*0.5; % [eV] energy scale of the problem
T0 = hbar/E0; % [s] time scale for the problem

n=1; % index of eigenstate
eigenstate_n_E=zeros(N,1);
eigenstate_n_E(n+1)=1;
pd_eigenstate_n_E=conj(eigenstate_n_E).*eigenstate_n_E; % possibility distribution
eigenstate_n_X=UXE*eigenstate_n_E;
pd_eigenstate_n_X=conj(eigenstate_n_X).*eigenstate_n_X;

t = linspace(0,TimeSpan*T0,nt); % [s] time vector

psi_0_case = 1;
psi_0 = zeros(N,1); % zero vector for initial state
switch psi_0_case
    case 1
        psi_0(1) = 1; % initial state is the ground state
    case 2
        psi_0(2) = 2; % initial state is the ground state
end

psi_t_E = zeros(N,nt); % time-varying state psi(t) [energy basis]
psi_t_X = zeros(N,nt); % time-varying state psi(t) [site basis]
for t_idx = 1:nt
    Ut = expm(-1i * H * t(t_idx)/hbar);
    psi_t_E(:,t_idx) = Ut * psi_0;
    psi_t_X(:,t_idx) = UXE * psi_t_E(:,t_idx);
end

% probability densities
pd_t_E = conj(psi_t_E).*psi_t_E;
pd_t_X = conj(psi_t_X).*psi_t_X;
%% Visualization
% Energy basis visualization
subplot(1,2,1);
bar3(pd_t_E);
grid on;
set(gca,'Fontsize',16,'Fontname','Times');
% xlabel('Probability','Interpreter','latex');
% ylabel('$n$','Interpreter','latex');

% site basis visualization
% subplot(1,3,2:3);
% bar(diag(Xeigval),pd_eigenstate_n_X);
% grid on;
% ylabel('Probability','Interpreter','latex');
% xlabel('$X(nm)$','Interpreter','latex');