%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%一个矩阵的求解
%备注：具有自旋轨道耦合的单粒子哈密顿量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%矩阵元涉及的变量定义
%%%%%%%%%%%%%%%%%%%%
%矩阵阶数：2N
N=500;
%自变量最大取值
rhomax=100;
%自变量最小取值
rhomin=0;
%difference value:DF
DF=(rhomax-rhomin)/(N-1);
%%%%%%%%%%%%%%%%%%
%自变量:rho=r的平方
rho=zeros(1,N);
%赋值
rho(1)=rhomin;
for i=2:N
    rho(i)=rho(i-1)+DF;
end
%%%%%%%%%%%%%%%%%%
%光强强度调节量：I0
I10=1;
I20=1;
%光场频率:w
w=5;
%拉盖尔高斯函数:L
L1=0.5;
L2=0.5;
%拉盖尔高斯主量子数：n
n1=1;
n2=-1;
%光场函数：I
I1=zeros(1,N);
I2=zeros(1,N);
for i=1:N
    I1(i)=I10*(rho(i)/w^2)^abs(n1)*(L1*exp(-rho(i)/w^2))^2;
    I2(i)=I20*(rho(i)/w^2)^abs(n2)*(L2*exp(-rho(i)/w^2))^2;
end
%拉曼耦合强度调节量：Omiga0
Omiga0=0;
%拉曼耦合项
Omiga=zeros(1,N);
for i=1:N
    Omiga(i)=Omiga0*sqrt((I10*(L1*exp(-rho(i)/w^2))^2)*(I20*(L2*exp(-rho(i)/w^2))^2));
end
%光场引入的等效失谐系数：X11,X12,X21,X22
X11=-1;
X12=-1;
X21=-1;
X22=-1;
%旋转波近似引入的失谐量：Det_RWA
Det_RWA=0;
%能级总失谐量：Det
Det=zeros(2,N);
for i=1:N
    Det(1,i)=X11*I1(i)+X21*I2(i)+Det_RWA;
    Det(2,i)=X12*I1(i)+X22*I2(i)-Det_RWA;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
%自变量：rho
%自变量：Lz
%%%%%%%%%%%%
%循环计数器：n
n=0;
%自变量：Lz=z方向的角动量
Lzmin=-5; 
Lzmax=5;
figure
hold on;
H=zeros(N*2,N*2);
for Lz=Lzmin:Lzmax
    %实验坐标系下的自旋角动量
    s1=abs(Lz-n1);
    s2=abs(Lz-n2);
    %为了方便标识s以及后面的计算
    s=zeros(2,N);
    for i=1:N
        s(1,i)=s1;
        s(2,i)=s2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %矩阵元分五部分：四部分+其它为0部分
    %矩阵元为作外积，因此j,k=|j,k>，jj,kk=<jj,kk|
    for k=1:2
        for kk=1:2
            for j=1:N
                for jj=1:N
                    
                    if j==jj&&k==kk
                        H(j+(k-1)*N,jj+(kk-1)*N)=4*rho(jj)/DF^2+2*(s(kk,jj)+1)/DF+Det(kk,jj)+rho(jj)/2;
                    end
                    if j==jj&&k==kk-1
                        H(j+(k-1)*N,jj+(kk-1)*N)=Omiga(jj)*rho(jj)^((abs(n1)+abs(n2)+s2-s1)/2);
                    end
                    if j==jj&&k==kk+1
                        H(j+(k-1)*N,jj+(kk-1)*N)=Omiga(jj)*rho(jj)^((abs(n1)+abs(n2)+s1-s2)/2);
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

    plot(Lz,EnergySpectrum(n,:),'s','MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',6);
end
xlabel('I=0');
xlim([-5,5]);
ylim([0,10]);
    

