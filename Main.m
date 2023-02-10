% This code was developed by Aref Miri Rekavandi for the paper: Rekavandi, A. M., Seghouane, A. K., & Evans, R. J. (2021). 
% Robust subspace detectors based on ?-divergence with application to detection in imaging. IEEE Transactions on Image Processing, 
% 30, 5017-5031. 
% If you use this code in your study, kindly cite the aforementioned paper.

clc
clear all
% close all
%% initialization for non-Gaussian case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indim=200;       %Dimension of input signal y
N=8000;          %Number of test samples
ttt=2;           % Number of colums in matrix H
ppp=2;           % Number of colums in matrix B
x=1:indim;
H=zeros(indim,ttt);
for i=1:ttt
  H(:,i)=cos(2*i*3.14*(1/indim)*x)';
end
H=H-mean(H);

B=zeros(indim,ppp);
for i=1:ppp
  B(:,i)=cos(2*(i-1)*3.14*(7/indim)*x)';
end

C=[H B];
mu1=[zeros(1,N/2) ones(1,N/2)];
alpha=0.75;
ASNR=-10;       % SNR
taw=1-alpha;
sw=1;           % sw=1 makes an asymmetric distribution
mue1=0;         % Noise mean in first component
mue2=2;         % Noise mean in second component
sigma1=1;       % Noise SD in first component
sigma2=4;       % Noise SD in second component
ep=0.15;        % Outlier rate
noise=noise_model(N,indim,ep,sigma1,mue1,sigma2,mue2,sw);
figure
[AAA BBB]=hist(noise(:),1000);
plot(BBB,AAA/N,'LineWidth',3)
xlabel('Noise value', 'Interpreter', 'LaTeX')
ylabel('Density', 'Interpreter', 'LaTeX')

 TETA1=random('uniform',0.3,0.300001,ttt,1);

 scale=sqrt(((indim*sigma1^2*10^(ASNR/10))/((H*TETA1)'*(H*TETA1))));
 SNR=10*log10((((scale*H*TETA1)'*(scale*H*TETA1))/indim)/(sigma1^2));
 
%% Making observations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:N
    TETA1=random('uniform',0.3,0.3001,ttt,1);
    X(:,k)=mu1(k)*H*TETA1;
    phi=random('uniform',0.3,0.3001,ppp,1);
    I(:,k)=B*phi;
    Y(:,k)=scale*X(:,k)+I(:,k)+noise(:,k);
end

%% Diffrent Tests with plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Tmsd TmsdAlpha Twald TraoAlpha THuber TwaldAlpha]=tests(H,C,indim,ttt,ppp,N,alpha,taw,Y);
figure
subplot(3,2,1)
 plot([Tmsd(1:N/2) zeros(1,N/2)])
 title('MSD')
 hold on
 plot([zeros(1,N/2) Tmsd(N/2+1:N)])
subplot(3,2,2)
 plot([TmsdAlpha(1:N/2) zeros(1,N/2)])
 title('Robust LRT')
 hold on
 plot([zeros(1,N/2) TmsdAlpha(N/2+1:N)])
subplot(3,2,3)
 plot([Twald(1:N/2) zeros(1,N/2)])
 title('Wald')
 hold on
 plot([zeros(1,N/2) Twald(N/2+1:N)])
subplot(3,2,4)
 plot([TwaldAlpha(1:N/2) zeros(1,N/2)])
 title('Robust Wald')
 hold on
 plot([zeros(1,N/2) TwaldAlpha(N/2+1:N)])
subplot(3,2,5)
 plot([THuber(1:N/2) zeros(1,N/2)])
 title('Rao Huber')
 hold on
 plot([zeros(1,N/2) THuber(N/2+1:N)])
subplot(3,2,6)
 plot([TraoAlpha(1:N/2) zeros(1,N/2)])
 title('Robust Rao')
 hold on
 plot([zeros(1,N/2) TraoAlpha(N/2+1:N)])

figure
h1 = histogram(TmsdAlpha(1:N/2),80);
hold on
h2 = histogram(TmsdAlpha(N/2+1:N),80);
legend({['$\mathcal{H}_0$'],['$\mathcal{H}_1$']}, ...
        'Interpreter', 'LaTeX')
    xlabel(['$\alpha$-LRT value'], 'Interpreter', 'LaTeX')
    ylabel('\#Sample', 'Interpreter', 'LaTeX')

%% ROC calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
last=max([max(Tmsd),max(TmsdAlpha), max(Twald), max(TwaldAlpha), max(THuber), max(TraoAlpha)]);
starting=min([min(Tmsd),min(TmsdAlpha), min(Twald), min(TwaldAlpha), min(THuber), min(TraoAlpha)]);
step=(last-starting)/1000000;
for th=starting:step:last
    pd1(i)=(200*(sum(Tmsd(N/2+1:N)>th)))/N;
    pf1(i)=(200*(sum(Tmsd(1:N/2)>th)))/N;
    
    pd2(i)=(200*(sum(TmsdAlpha(N/2+1:N)>th)))/N;
    pf2(i)=(200*(sum(TmsdAlpha(1:N/2)>th)))/N;  
    
    pd3(i)=(200*(sum(Twald(N/2+1:N)>th)))/N;
    pf3(i)=(200*(sum(Twald(1:N/2)>th)))/N; 
    
    pd4(i)=(200*(sum(TwaldAlpha(N/2+1:N)>th)))/N;
    pf4(i)=(200*(sum(TwaldAlpha(1:N/2)>th)))/N; 
    
    pd5(i)=(200*(sum(THuber(N/2+1:N)>th)))/N;
    pf5(i)=(200*(sum(THuber(1:N/2)>th)))/N; 
    
    pd6(i)=(200*(sum(TraoAlpha(N/2+1:N)>th)))/N;
    pf6(i)=(200*(sum(TraoAlpha(1:N/2)>th)))/N; 
    i=i+1;
end

%% initialization for Gaussian case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sw=0;
mue1=0;
mue2=5;
sigma1=1;
sigma2=0.1;
ep=0;
noise=noise_model(N,indim,ep,sigma1,mue1,sigma2,mue2,sw);
%% Making observations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:N
    TETA1=random('uniform',0.3,0.3001,ttt,1);
    X(:,k)=mu1(k)*H*TETA1;
    phi=random('uniform',0.3,0.3001,ppp,1);
    I(:,k)=B*phi;
    Y(:,k)=scale*X(:,k)+I(:,k)+noise(:,k);
end
%% Diffrent Tess %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Tmsd1 TmsdAlpha1 Twald1 TraoAlpha1 THuber1 TwaldAlpha1]=tests(H,C,indim,ttt,ppp,N,alpha,taw,Y);

%% ROC calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
last=max([max(Tmsd1),max(TmsdAlpha1), max(Twald1), max(TwaldAlpha1), max(THuber1), max(TraoAlpha1)]);
starting=min([min(Tmsd1),min(TmsdAlpha1), min(Twald1), min(TwaldAlpha1), min(THuber1), min(TraoAlpha1)]);
step=(last-starting)/1000000;
for th=starting:step:last
    pd11(i)=(200*(sum(Tmsd1(N/2+1:N)>th)))/N;
    pf11(i)=(200*(sum(Tmsd1(1:N/2)>th)))/N;
      
    pd22(i)=(200*(sum(TmsdAlpha1(N/2+1:N)>th)))/N;
    pf22(i)=(200*(sum(TmsdAlpha1(1:N/2)>th)))/N; 
    
    pd33(i)=(200*(sum(Twald1(N/2+1:N)>th)))/N;
    pf33(i)=(200*(sum(Twald1(1:N/2)>th)))/N;
    
    pd44(i)=(200*(sum(TwaldAlpha1(N/2+1:N)>th)))/N;
    pf44(i)=(200*(sum(TwaldAlpha1(1:N/2)>th)))/N; 
    
    pd55(i)=(200*(sum(THuber1(N/2+1:N)>th)))/N;
    pf55(i)=(200*(sum(THuber1(1:N/2)>th)))/N; 
    
    pd66(i)=(200*(sum(TraoAlpha1(N/2+1:N)>th)))/N;
    pf66(i)=(200*(sum(TraoAlpha1(1:N/2)>th)))/N; 
    i=i+1;
end
%% ROC plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(pf11,pd11,'b','LineWidth',2)
hold on
plot(pf22,pd22,'r','LineWidth',2)
hold on
plot(pf33,pd33,'g','LineWidth',2)
hold on
plot(pf44,pd44,'k','LineWidth',2)
hold on
plot(pf55,pd55,'c','LineWidth',2)
hold on
plot(pf66,pd66,'m','LineWidth',2)
hold on
plot(pf1,pd1,'--b','LineWidth',2)
hold on
plot(pf2,pd2,'--r','LineWidth',2)
hold on
plot(pf3,pd3,'--g','LineWidth',2)
hold on
plot(pf4,pd4,'--k','LineWidth',2)
hold on
plot(pf5,pd5,'--c','LineWidth',2)
hold on
plot(pf6,pd6,'--m','LineWidth',2)

grid on

    legend({'MSD in Gaussian noise',['Robust LRT in Gaussian noise with $\alpha$=' num2str(alpha) ''],'Wald in Gaussian noise',['Robust Wald in Gaussian noise with $\alpha$=' num2str(alpha) ''],'Huber based Rao in  Gaussian noise',['Robust Rao in Gaussian noise with $\alpha$=' num2str(alpha) ''],'MSD in contaminated Gaussian noise',['Robust LRT in contaminated Gaussian noise with $\alpha$=' num2str(alpha) ''],'Wald in contaminated Gaussian noise',['Robust Wald in contaminated Gaussian noise with $\alpha$=' num2str(alpha) ''],'Huber based Rao in contaminated Gaussian noise',['Robust Rao in contaminated Gaussian noise with $\alpha$=' num2str(alpha) '']}, ...
        'Interpreter', 'LaTeX')
    xlabel('Probability of False Alarm (\%)', 'Interpreter', 'LaTeX')
    ylabel('Probability of Detection (\%)', 'Interpreter', 'LaTeX')
    title('(b)', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')

