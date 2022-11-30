% ====================== 单粒子 SOAMC =================
% ====================== 一个矩阵的求解 ================
% ====================== 备注：具有自旋轨道耦合的单粒子哈密顿量 ===========

% ====================== 矩阵元涉及的变量定义 ============
close all; clear all;
N=500; % 矩阵阶数：2N
rhomax=100; % 自变量最大取值
rhomin=0; % 自变量最小取值

% === 自变量：rho = r^2，柱坐标下 r = (r) ============
rho = linspace(rhomin, rhomax, N); % 这样写效果不一样吗
DF = rho(2) - rho(1); 

% ============ 光强强度调节量：I0 =======================
I0 = 0;
I10=I0;
I20=I0;
w=5; % 光场频率：w，束宽

L1=0.5; % 拉盖尔高斯函数：L，常数
L2=0.5;

n1=1; % 拉盖尔高斯主量子数：n，也就是 phase winding
n2=-1;
%光场函数：I
I1 = I10*(rho/w^2).^abs(n1).*(L1*exp(-rho/w^2)).^2;
I2 = I20*(rho/w^2).^abs(n2).*(L2*exp(-rho/w^2)).^2;

omega0 = 5; % 拉曼耦合强度调节量：Omega0
omega = omega0*sqrt(I1.*I2);
% ================== 光场引入的等效失谐系数：X11,X12,X21,X22 ===============
X11 = -1; X12 = -1;
X21 = -1; X22 = -1;

Det_RWA = 0; % 旋转波近似引入的失谐量：Det_RWA
Det = zeros(2,N); % 能级总失谐量：Det
Det(2,:) = X12*I1+X22*I2 - Det_RWA;
Det(1,:) = X11*I1+X21*I2 + Det_RWA;

% ===================== 矩阵 =======================
% ===================== 自变量：rho ================
% ===================== 自变量：Lz =================
% ===================== 循环计数器：n ===============
n = 0;
% ===================== 自变量：Lz = z (z 方向的角动量) ==========
Lzmin = -5;  Lzmax = 5;
figure
hold on;
H = zeros(N*2, N*2);
for Lz = Lzmin : Lzmax
    s1 = abs(Lz - n1); % 实验坐标系下的自旋角动量
    s2 = abs(Lz - n2);
    s = zeros(2, N); % 为了方便标识 s 以及后面的计算
    for i=1:N
        s(1,i)=s1;
        s(2,i)=s2;
    end
    % =========================================
    % ================= 矩阵元分五部分：四部分+其它为0部分 ===========
    % ================= 矩阵元为作外积，因此j,k=|j,k>，jj,kk=<jj,kk|
    for k=1:2
        for kk=1:2
            for j=1:N
                for jj=1:N
                    
                    if j==jj&&k==kk
                        H(j+(k-1)*N,jj+(kk-1)*N)=4*rho(jj)/DF^2+2*(s(kk,jj)+1)/DF+Det(kk,jj)+rho(jj)/2;
                    end
                    if j==jj&&k==kk-1
                        H(j+(k-1)*N,jj+(kk-1)*N)=omega(jj)*rho(jj)^((abs(n1)+abs(n2)+s2-s1)/2);
                    end
                    if j==jj&&k==kk+1
                        H(j+(k-1)*N,jj+(kk-1)*N)=omega(jj)*rho(jj)^((abs(n1)+abs(n2)+s1-s2)/2);
                    end
                    
                    if j~=1&&j~=N
                        if j==jj+1&&k==kk
                            H(j+(k-1)*N,jj+(kk-1)*N)=-2*rho(jj)/DF^2-2*(s(kk,jj)+1)/DF;
                        end
                        if j==jj-1&&k==kk
                            H(j+(k-1)*N,jj+(kk-1)*N)=-2*rho(jj)/DF^2;
                        end
                    end
                    
                    if j==1
                        %因为jj~=0
                        %if j==jj+1&&k==kk
                        %    H(j+(k-1)*N,jj+(kk-1)*N)=-2*rho(jj)/DF^2-2*(s(kk,jj)+1)/DF;
                        %end
                        if j==jj-1&&k==kk
                            H(j+(k-1)*N,jj+(kk-1)*N)=-2*rho(jj)/DF^2;
                        end
                    end
                    if j==N
                        if j==jj+1&&k==kk
                            H(j+(k-1)*N,jj+(kk-1)*N)=-2*rho(jj)/DF^2-2*(s(kk,jj)+1)/DF;
                        end
                        %因为jj~=N+1
                        %if j==jj-1&&k==kk
                        %    H(j+(k-1)*N,jj+(kk-1)*N)=-2*rho(jj)/DF^2;
                        %end
                    end
                end
            end
        end
    end
    n=n+1;
    [Vec,Eig]=eig(H);
    [D_sort,index] = sort(diag(Eig));
    EnergySpectrum(n,:)=D_sort(1:15);

    plot(Lz,EnergySpectrum(n,:),'s','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',5);
end
axis([-5 5 0 20]);
xlabel('\it L_z', 'FontName','Times New Roman');
ylabel('\it E', 'FontName','Times New Roman');

