clear;clc;
for u=1:1
if(u==1)k=['rzdata','.out'];end
% if(u==2)k=['rzdata_ion_Er','.out'];end
% if(u==3)k=['rzdata_ele_noEr','.out'];end
% if(u==4)k=['rzdata_ele_Er','.out'];end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%                  read data                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    atmp=load(k);
    mn=size(atmp);m=mn(1);n=mn(2);
    mpsi=n;mtor=32;
%     mtime=m/9;
    rhoi=1;%2.13e-4/2.955e-4;
%     R0=5617.9; %major radius unit mm
    timestart=0;timeend=0;%mtime-1;
     N=8;% N is range of Cr grid
     Cr_store_timeave=zeros(N,1);
    for S=timestart:timeend
    K=mtor*S;
    a=atmp(1+K:mtor+K,:);
    Cr_store=zeros(N,1);
    for k=1:mtor
    a1=a(k,:);
    psi0=20;psi1=mpsi-20;
    Mpsi=psi1-psi0+1;
   
    Cr=zeros(N,1);
    Cr1=zeros(N,1);Cr2=zeros(N,1);
    for j=1:N
        jj=j-1;
        for i=psi0:1:psi1
            ii=i-psi0+1;
           % p1=(a1(i-jj)+a1(i+jj))/2;
            p1=a1(i-jj);
            Cr(j)=Cr(j)+a1(i)*p1;
            Cr1(j)=Cr1(j)+a1(i)*a1(i);
            Cr2(j)=Cr2(j)+p1^2;
        end
    end
    Cr=Cr./sqrt(Cr1.*Cr2);
    Cr_store=Cr_store+Cr;
%      plot(0:N-1,Cr,'LineWidth',2);hold on;
    end
    Cr_store=Cr_store/mtor;
    x=(0:N-1)';x=x/rhoi;
%    if(u==1) plot(x,Cr_store,'b-','LineWidth',2);hold on;end
%    if(u==2) plot(x,Cr_store,'r-','LineWidth',2);hold on;end
    
%     mean(x.*Cr_store)/mean(Cr_store)

   Cr_store_timeave=Cr_store_timeave+Cr_store/(timeend-timestart+1);

% surf(1:mpsi,1:mtor,a);
% axis tight;view(0,90) ;
% shading interp;
% xlabel('radial grid');ylabel('parallel grid');
    end
   plot(x,Cr_store_timeave,'-','LineWidth',2);hold on;
   
% mean(x.*Cr_store_timeave)/mean(Cr_store_timeave)
end
% mean(x.*Cr_store_timeave)/mean(Cr_store_timeave)
% plot([0,12],[0.3679,0.3679],'k--');
%fit exp initial drop
%    Lr=1.7;
%  

% Lr=2.6;
%    r=0:0.1:N-1;r=r/rhoi;
%    fity=exp(-r/Lr);
%    plot(r,fity,'b--','LineWidth',2);hold on;
% Lr=1.61;
%    r=0:N-1;r=r/rhoi;
%    fity=exp(-r/Lr);
%    plot(r,fity,'r--','LineWidth',2);hold on;
set(gca,'FontSize',24);
grid on;
% plot([0,12],[0.3679,0.3679],'k--');
% legend('Cr(\Deltar) noE_e_q','Cr(\Deltar)     E_e_q');
%  legend('noE_e_q','    E_e_q');
% legend('boxoff');
xlabel('\Deltar(\rho_i)');ylabel('C_r');
% % % xlabel('radial grid');ylabel('Corelation length(grid)');
% % % set(gca,'YLim',[0 12]);
% set(gca,'XLim',[0 6]);
set(gca,'YLim',[0 1]);
set (gcf,'Position',[100,50,700,600], 'color','w');




