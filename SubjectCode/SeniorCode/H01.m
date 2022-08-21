%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%һ����������
%��ע���������������ϵĵ����ӹ��ܶ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%����Ԫ�漰�ı�������
%%%%%%%%%%%%%%%%%%%%
%���������2N
N=500;
%�Ա������ȡֵ
rhomax=100;
%�Ա�����Сȡֵ
rhomin=0;
%difference value:DF
DF=(rhomax-rhomin)/(N-1);
%%%%%%%%%%%%%%%%%%
%�Ա���:rho=r��ƽ��
rho=zeros(1,N);
%��ֵ
rho(1)=rhomin;
for i=2:N
    rho(i)=rho(i-1)+DF;
end
%%%%%%%%%%%%%%%%%%
%��ǿǿ�ȵ�������I0
I10=1;
I20=1;
%�ⳡƵ��:w
w=5;
%���Ƕ���˹����:L
L1=0.5;
L2=0.5;
%���Ƕ���˹����������n
n1=1;
n2=-1;
%�ⳡ������I
I1=zeros(1,N);
I2=zeros(1,N);
for i=1:N
    I1(i)=I10*(rho(i)/w^2)^abs(n1)*(L1*exp(-rho(i)/w^2))^2;
    I2(i)=I20*(rho(i)/w^2)^abs(n2)*(L2*exp(-rho(i)/w^2))^2;
end
%�������ǿ�ȵ�������Omiga0
Omiga0=0;
%���������
Omiga=zeros(1,N);
for i=1:N
    Omiga(i)=Omiga0*sqrt((I10*(L1*exp(-rho(i)/w^2))^2)*(I20*(L2*exp(-rho(i)/w^2))^2));
end
%�ⳡ����ĵ�Чʧгϵ����X11,X12,X21,X22
X11=-1;
X12=-1;
X21=-1;
X22=-1;
%��ת�����������ʧг����Det_RWA
Det_RWA=0;
%�ܼ���ʧг����Det
Det=zeros(2,N);
for i=1:N
    Det(1,i)=X11*I1(i)+X21*I2(i)+Det_RWA;
    Det(2,i)=X12*I1(i)+X22*I2(i)-Det_RWA;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
%�Ա�����rho
%�Ա�����Lz
%%%%%%%%%%%%
%ѭ����������n
n=0;
%�Ա�����Lz=z����ĽǶ���
Lzmin=-5; 
Lzmax=5;
figure
hold on;
H=zeros(N*2,N*2);
for Lz=Lzmin:Lzmax
    %ʵ������ϵ�µ������Ƕ���
    s1=abs(Lz-n1);
    s2=abs(Lz-n2);
    %Ϊ�˷����ʶs�Լ�����ļ���
    s=zeros(2,N);
    for i=1:N
        s(1,i)=s1;
        s(2,i)=s2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %����Ԫ���岿�֣��Ĳ���+����Ϊ0����
    %����ԪΪ����������j,k=|j,k>��jj,kk=<jj,kk|
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
                        %��Ϊjj~=0
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
                        %��Ϊjj~=N+1
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
    

