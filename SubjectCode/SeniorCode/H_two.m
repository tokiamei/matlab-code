%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spin-angular momentum coupling in species A and non spin-angular momentum
%coupling in species B; through spin-exchange interaction, the spin and
%angular momentum degree of freedom in species B are effectively coupled.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%two D imaginary time evolution by which the ground state is obtained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% program starts here
%%parameterize the discrete grids in the field
spalength=5;
ax=-spalength;bx=spalength;
ay=-spalength;by=spalength;
N=100;
Nx=N;Ny=N;
x=linspace(ax,bx,Nx);y=linspace(ay,by,Ny);
h=x(2)-x(1);
Maxerr=1e-6;  % the maximum error
itererr=1.0;  % convergence criterion parameter
dt=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%parameterize characteristics for species A and B

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NA£¬NB
% atom number
NA=1000;
NB=1000; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MA=1;MB=1; % atom mass
omegaA=1;omegaB=1; %trap frequency

g11=0.01e0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%g12
g12=0.005e0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g22=0.01e0; % nonlinear interaction parameter for species A

h11=0.01e0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%h12
h12=0.015e0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h22=0.01e0; % nonlinear interaction parameter for species B

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%parameterize the interaction parameters between species A and B
gamma=0.01e0; % density-density interaction and

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%beta
beta=0.02e0;  %spin-exchange interaction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%parameterize the LG laser
w=2.0; % the width of waist of the laser
l1=1; % the angular momentum quantum number of LG1
l2=-1; % the angular momentum quantum number of LG2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OmegaR
OmegaR=2; % raman frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta=0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%intermediate variables
RNMAB=NB*MA/(NA*MB);
ROmegasquare=NB*MB*omegaB^2/(NA*MA*omegaA^2);
[xx,yy]=meshgrid(x,y);
Omega1=OmegaR*exp(-1i*(l1-l2)*atan2(yy,xx)).*(sqrt(xx.^2+yy.^2)./w).^(abs(l1)+abs(l2)).*exp(-2*(xx.^2+yy.^2)./w^2);
Omega2=OmegaR*exp(1i*(l1-l2)*atan2(yy,xx)).*(sqrt(xx.^2+yy.^2)./w).^(abs(l1)+abs(l2)).*exp(-2*(xx.^2+yy.^2)./w^2);

chi11=NA*g11;chi12=NA*g12;chi22=NA*g22;
kapa11=NB^2/NA*h11;kapa12=NB^2/NA*h12;kapa22=NB^2/NA*h22;
Idd=NB*gamma;Iss=NB*beta;

COA=1./(1+2*dt/h^2+dt*(xx.^2+yy.^2)/2);
COB=1./(1+2*dt*RNMAB/h^2+dt*ROmegasquare*(xx.^2+yy.^2)/2);
COAK=dt/(2*h^2);COBK=dt*RNMAB/(2*h^2);
CSOA1=dt*Omega1;CSOA2=dt*Omega2;
Cdelta=dt*delta;
CA11=dt*chi11;CA12=dt*chi12;CA22=dt*chi22;
CB11=dt*kapa11;CB12=dt*kapa12;CB22=dt*kapa22;
Cdd=dt*Idd;Css=dt*Iss;

ECA=1/(2*h^2);
TrapA1=(xx.^2+yy.^2)/2+delta; TrapA2=(xx.^2+yy.^2)/2-delta;
ECB=RNMAB/(2*h^2);
TrapB=ROmegasquare*(xx.^2+yy.^2)/2;
EA11=chi11/2;EA12=chi12;EA22=chi22/2;
EB11=kapa11/2;EB12=kapa12;EB22=kapa22/2;

clear xx yy;

%%% initialize the wavefunction
psiAup = rand(Nx,Ny)+1i*rand(Nx,Ny);
psiAdn = rand(Nx,Ny)+1i*rand(Nx,Ny);
psiBup = rand(Nx,Ny)+1i*rand(Nx,Ny);
psiBdn = rand(Nx,Ny)+1i*rand(Nx,Ny);
psiAup(1,:)=0;psiAup(Nx,:)=0;
psiAup(:,1)=0;psiAup(:,Ny)=0;
psiAdn(1,:)=0;psiAdn(Nx,:)=0;
psiAdn(:,1)=0;psiAdn(:,Ny)=0;
psiBup(1,:)=0;psiBup(Nx,:)=0;
psiBup(:,1)=0;psiBup(:,Ny)=0;
psiBdn(1,:)=0;psiBdn(Nx,:)=0;
psiBdn(:,1)=0;psiBdn(:,Ny)=0;
% the output wf during iteration
psiAuptemp2=zeros(Nx,Ny);
psiAdntemp2=zeros(Nx,Ny);

psiBuptemp2=zeros(Nx,Ny);
psiBdntemp2=zeros(Nx,Ny);

% the input wf during iteration
psiAuptemp1=zeros(Nx,Ny);
psiAdntemp1=zeros(Nx,Ny);

psiBuptemp1=zeros(Nx,Ny);
psiBdntemp1=zeros(Nx,Ny);

%%% normalize the wavefunction
NormA=h^2*(sum(abs(psiAup).^2+abs(psiAdn).^2,'all'));
NormB=h^2*(sum(abs(psiBup).^2+abs(psiBdn).^2,'all'));
psiAup = psiAup/sqrt(NormA);
psiAdn = psiAdn/sqrt(NormA);
psiBup = psiBup/sqrt(NormB);
psiBdn = psiBdn/sqrt(NormB);

T=0;TK=0;
energy=zeros(1,1000000);

while (itererr>Maxerr)
    TK=TK+1;T=dt*TK
%    psiAuptemp2=solveAup(Nx,Ny,Maxerr,COA,COAK,CSOA1,Cdelta,CA11,CA12,Cdd,Css,psiAup,psiAdn,psiBup,psiBdn,psiAuptemp1);
    diff=1.d0;
    psiAuptemp1=psiAup;
    psiAuptemp2=psiAuptemp1;
      while (diff>Maxerr)
         psiAuptemp2(2:Nx-1,2:Ny-1)=COA(2:Nx-1,2:Ny-1).*(psiAup(2:Nx-1,2:Ny-1)+COAK*(psiAuptemp1(3:Nx,2:Ny-1)+psiAuptemp1(1:Nx-2,2:Ny-1)+psiAuptemp1(2:Nx-1,3:Ny)+psiAuptemp1(2:Nx-1,1:Ny-2))...
         -CSOA1(2:Nx-1,2:Ny-1).*psiAdn(2:Nx-1,2:Ny-1)-Cdelta*psiAuptemp1(2:Nx-1,2:Ny-1)-CA11*abs(psiAup(2:Nx-1,2:Ny-1)).^2.*psiAuptemp1(2:Nx-1,2:Ny-1)-CA12*abs(psiAdn(2:Nx-1,2:Ny-1)).^2.*psiAuptemp1(2:Nx-1,2:Ny-1)...
         -Cdd*(abs(psiBup(2:Nx-1,2:Ny-1)).^2+abs(psiBdn(2:Nx-1,2:Ny-1)).^2).*psiAuptemp1(2:Nx-1,2:Ny-1)-Css*conj(psiBdn(2:Nx-1,2:Ny-1)).*psiBup(2:Nx-1,2:Ny-1).*psiAdn(2:Nx-1,2:Ny-1));
         diff=max(max(abs(psiAuptemp2-psiAuptemp1)));
         psiAuptemp1=psiAuptemp2;
       end
%    psiAdntemp2=solveAdn(Nx,Ny,Maxerr,COA,COAK,CSOA2,Cdelta,CA22,CA12,Cdd,Css,psiAup,psiAdn,psiBup,psiBdn,psiAdntemp1);
    diff=1.d0;
    psiAdntemp1=psiAdn;
    psiAdntemp2=psiAdntemp1;
    while (diff>Maxerr)
      psiAdntemp2(2:Nx-1,2:Ny-1)=COA(2:Nx-1,2:Ny-1).*(psiAdn(2:Nx-1,2:Ny-1)+COAK*(psiAdntemp1(3:Nx,2:Ny-1)+psiAdntemp1(1:Nx-2,2:Ny-1)+psiAdntemp1(2:Nx-1,3:Ny)+psiAdntemp1(2:Nx-1,1:Ny-2))...
      -CSOA2(2:Nx-1,2:Ny-1).*psiAup(2:Nx-1,2:Ny-1)+Cdelta*psiAdntemp1(2:Nx-1,2:Ny-1)-CA22*abs(psiAdn(2:Nx-1,2:Ny-1)).^2.*psiAdntemp1(2:Nx-1,2:Ny-1)-CA12*abs(psiAup(2:Nx-1,2:Ny-1)).^2.*psiAdntemp1(2:Nx-1,2:Ny-1)...
      -Cdd*(abs(psiBup(2:Nx-1,2:Ny-1)).^2+abs(psiBdn(2:Nx-1,2:Ny-1)).^2).*psiAdntemp1(2:Nx-1,2:Ny-1)-Css*conj(psiBup(2:Nx-1,2:Ny-1)).*psiBdn(2:Nx-1,2:Ny-1).*psiAup(2:Nx-1,2:Ny-1));
      diff=max(max(abs(psiAdntemp2-psiAdntemp1)));
      psiAdntemp1=psiAdntemp2;
    end
%    psiBuptemp2=solveBup(Nx,Ny,Maxerr,COB,COBK,CB11,CB12,Cdd,Css,psiAup,psiAdn,psiBup,psiBdn,psiBuptemp1);
    diff=1.d0;
    psiBuptemp1=psiBup;
    psiBuptemp2=psiBuptemp1;
    while (diff>Maxerr)
      psiBuptemp2(2:Nx-1,2:Ny-1)=COB(2:Nx-1,2:Ny-1).*(psiBup(2:Nx-1,2:Ny-1)+COBK*(psiBuptemp1(3:Nx,2:Ny-1)+psiBuptemp1(1:Nx-2,2:Ny-1)+psiBuptemp1(2:Nx-1,3:Ny)+psiBuptemp1(2:Nx-1,1:Ny-2))...
      -CB11*abs(psiBup(2:Nx-1,2:Ny-1)).^2.*psiBuptemp1(2:Nx-1,2:Ny-1)-CB12*abs(psiBdn(2:Nx-1,2:Ny-1)).^2.*psiBuptemp1(2:Nx-1,2:Ny-1)...
      -Cdd*(abs(psiAup(2:Nx-1,2:Ny-1)).^2+abs(psiAdn(2:Nx-1,2:Ny-1)).^2).*psiBuptemp1(2:Nx-1,2:Ny-1)-Css*conj(psiAdn(2:Nx-1,2:Ny-1)).*psiAup(2:Nx-1,2:Ny-1).*psiBdn(2:Nx-1,2:Ny-1));
      diff=max(max(abs(psiBuptemp2-psiBuptemp1)));
      psiBuptemp1=psiBuptemp2;
    end
%    psiBdntemp2=solveBdn(Nx,Ny,Maxerr,COB,COBK,CB22,CB12,Cdd,Css,psiAup,psiAdn,psiBup,psiBdn,psiBdntemp1);
    diff=1.d0;
    psiBdntemp1=psiBdn;
    psiBdntemp2=psiBdntemp1;
    while (diff>Maxerr)
      psiBdntemp2(2:Nx-1,2:Ny-1)=COB(2:Nx-1,2:Ny-1).*(psiBdn(2:Nx-1,2:Ny-1)+COBK*(psiBdntemp1(3:Nx,2:Ny-1)+psiBdntemp1(1:Nx-2,2:Ny-1)+psiBdntemp1(2:Nx-1,3:Ny)+psiBdntemp1(2:Nx-1,1:Ny-2))...
      -CB22*abs(psiBdn(2:Nx-1,2:Ny-1)).^2.*psiBdntemp1(2:Nx-1,2:Ny-1)-CB12*abs(psiBup(2:Nx-1,2:Ny-1)).^2.*psiBdntemp1(2:Nx-1,2:Ny-1)...
      -Cdd*(abs(psiAup(2:Nx-1,2:Ny-1)).^2+abs(psiAdn(2:Nx-1,2:Ny-1)).^2).*psiBdntemp1(2:Nx-1,2:Ny-1)-Css*conj(psiAup(2:Nx-1,2:Ny-1)).*psiAdn(2:Nx-1,2:Ny-1).*psiBup(2:Nx-1,2:Ny-1));
      diff=max(max(abs(psiBdntemp2-psiBdntemp1)));
      psiBdntemp1=psiBdntemp2;
    end
    
    NormA=h^2*(sum(abs(psiAuptemp2).^2+abs(psiAdntemp2).^2,'all'));
    NormB=h^2*(sum(abs(psiBuptemp2).^2+abs(psiBdntemp2).^2,'all'));
    psiAuptemp2 = psiAuptemp2/sqrt(NormA);
    psiAdntemp2 = psiAdntemp2/sqrt(NormA);
    psiBuptemp2 = psiBuptemp2/sqrt(NormB);
    psiBdntemp2 = psiBdntemp2/sqrt(NormB);
    
    itererr=max(max(abs([psiAuptemp2-psiAup,psiAdntemp2-psiAdn,psiBuptemp2-psiBup,psiBdntemp2-psiBdn])));
    
    psiAup=psiAuptemp2;
    psiAdn=psiAdntemp2;
    psiBup=psiBuptemp2;
    psiBdn=psiBdntemp2;
    
    %% calculate the total energy!
    energy(TK)=h^2*sum(sum(-ECA*conj(psiAup(2:Nx-1,2:Ny-1)).*(psiAup(3:Nx,2:Ny-1)+psiAup(1:Nx-2,2:Ny-1)+psiAup(2:Nx-1,3:Ny)+psiAup(2:Nx-1,1:Ny-2)-4*psiAup(2:Nx-1,2:Ny-1))+TrapA1(2:Nx-1,2:Ny-1).*abs(psiAup(2:Nx-1,2:Ny-1)).^2 ...
                       -ECA*conj(psiAdn(2:Nx-1,2:Ny-1)).*(psiAdn(3:Nx,2:Ny-1)+psiAdn(1:Nx-2,2:Ny-1)+psiAdn(2:Nx-1,3:Ny)+psiAdn(2:Nx-1,1:Ny-2)-4*psiAdn(2:Nx-1,2:Ny-1))+TrapA2(2:Nx-1,2:Ny-1).*abs(psiAdn(2:Nx-1,2:Ny-1)).^2 ...
                       +Omega2(2:Nx-1,2:Ny-1).*conj(psiAdn(2:Nx-1,2:Ny-1)).*psiAup(2:Nx-1,2:Ny-1)+Omega1(2:Nx-1,2:Ny-1).*conj(psiAup(2:Nx-1,2:Ny-1)).*psiAdn(2:Nx-1,2:Ny-1) ...
                       -ECB*conj(psiBup(2:Nx-1,2:Ny-1)).*(psiBup(3:Nx,2:Ny-1)+psiBup(1:Nx-2,2:Ny-1)+psiBup(2:Nx-1,3:Ny)+psiBup(2:Nx-1,1:Ny-2)-4*psiBup(2:Nx-1,2:Ny-1))+TrapB(2:Nx-1,2:Ny-1).*abs(psiBup(2:Nx-1,2:Ny-1)).^2 ...
                       -ECB*conj(psiBdn(2:Nx-1,2:Ny-1)).*(psiBdn(3:Nx,2:Ny-1)+psiBdn(1:Nx-2,2:Ny-1)+psiBdn(2:Nx-1,3:Ny)+psiBdn(2:Nx-1,1:Ny-2)-4*psiBdn(2:Nx-1,2:Ny-1))+TrapB(2:Nx-1,2:Ny-1).*abs(psiBdn(2:Nx-1,2:Ny-1)).^2 ...
                       +EA11*abs(psiAup(2:Nx-1,2:Ny-1)).^4+EA12*abs(psiAup(2:Nx-1,2:Ny-1)).^2.*abs(psiAdn(2:Nx-1,2:Ny-1)).^2+EA22*abs(psiAdn(2:Nx-1,2:Ny-1)).^4+EB11*abs(psiBup(2:Nx-1,2:Ny-1)).^4 ...
                       +EB12*abs(psiBup(2:Nx-1,2:Ny-1)).^2.*abs(psiBdn(2:Nx-1,2:Ny-1)).^2+EB22*abs(psiBdn(2:Nx-1,2:Ny-1)).^4 ...
                       +Idd*(abs(psiAup(2:Nx-1,2:Ny-1)).^2+abs(psiAdn(2:Nx-1,2:Ny-1)).^2).*(abs(psiBup(2:Nx-1,2:Ny-1)).^2+abs(psiBdn(2:Nx-1,2:Ny-1)).^2) ...
                       +Iss*(conj(psiAup(2:Nx-1,2:Ny-1)).*conj(psiBdn(2:Nx-1,2:Ny-1)).*psiBup(2:Nx-1,2:Ny-1).*psiAdn(2:Nx-1,2:Ny-1)+conj(psiAdn(2:Nx-1,2:Ny-1)).*conj(psiBup(2:Nx-1,2:Ny-1)).*psiBdn(2:Nx-1,2:Ny-1).*psiAup(2:Nx-1,2:Ny-1))));


    
end

figure
subplot(3,4,1);surf(x,y,abs(psiAup).^2);view(2);shading interp;
subplot(3,4,2);surf(x,y,abs(psiAdn).^2);view(2);shading interp;
subplot(3,4,5);surf(x,y,angle(psiAup));view(2);shading interp;
subplot(3,4,6);surf(x,y,angle(psiAdn));view(2);shading interp;

subplot(3,4,3);surf(x,y,abs(psiBup).^2);view(2);shading interp;
subplot(3,4,4);surf(x,y,abs(psiBdn).^2);view(2);shading interp;
subplot(3,4,7);surf(x,y,angle(psiBup));view(2);shading interp;
subplot(3,4,8);surf(x,y,angle(psiBdn));view(2);shading interp;  

subplot(3,4,[9 10 11 12]);plot(linspace(0,T,TK),real(energy(1:TK)));
   