
function [Tmsd TmsdAlpha Twald TraoAlpha THuber TwaldAlpha]=tests(H,Htilda,indim,ttt,ppp,N,alpha,taw,Y);

PH=H*((H'*H)^(-1))*H';
S=Htilda(:,ttt+1:ppp+ttt);
PSO=eye(indim)-(([S])*inv((([S])')*([S]))*(([S])'));
G=PSO*H;
PG=((G)*inv(((G)')*(G))*((G)'));
PGO=eye(indim)-((G)*inv(((G)')*(G))*((G)'));


for i=1:N
    clc
    fprintf(' %d out of 100 is done!',round(100*i/N));

temp=Y(:,i);



T1=sum((PG*PSO*temp).^2);
T2=sum((PGO*PSO*temp).^2);

Tmsd(i)=T1/T2;

%%%%%%   Version1  %%%%%%%%%%
%%%%% Hypothesis 0 %%%%%%%%%%%
  W0=eye(indim);
  for l=1:30
      par0=((S'*W0*S)^(-1))*S'*W0*temp;
      w0=sum(W0);
      var0=(w0*((temp-S*par0).^2))/(sum(w0));
      for h=1:indim
          W0(h,h)=exp(-0.5*taw*((temp(h)-S(h,:)*par0)^2)/var0);
      end
  end 
  
 %%%%% Hypothesis 1 %%%%%%%%%%% 
  W1=eye(indim);
  for l=1:30
      par1=((Htilda'*W1*Htilda)^(-1))*Htilda'*W1*temp;
      w1=sum(W1);
      var1=(w1*((temp-Htilda*par1).^2))/(sum(w1));
      for h=1:indim
          W1(h,h)=exp(-0.5*taw*((temp(h)-Htilda(h,:)*par1)^2)/var1);
      end
  end 
 pltvar1(i)=var1;
 TmsdAlpha(i)=0; 
 for time=1:indim   
    TmsdAlpha(i)=TmsdAlpha(i)+((w1(time)/(var1^(taw/2)))-(w0(time)/(var0^(taw/2))));
 end
TmsdAlpha(i)=(2/(indim*taw*((2*pi)^(taw/2))))*TmsdAlpha(i);

%%%%%%% Hypothesis 1 Wald %%%%%%%%%%%%%%%%%
  par1=((Htilda'*Htilda)^(-1))*Htilda'*temp;
  var1=sum(((temp-Htilda*par1).^2))/indim;
  covinv=(Htilda'*Htilda);
  whitenedpar=covinv^(0.5)*par1;
  targetpar=whitenedpar(1:ttt);
  T1=(targetpar)'*(targetpar);
  T2=var1;
  Twald(i)=T1/T2;
  
%%%%% Hypothesis 0 Robust Rao %%%%%%%%%%%
  W0=eye(indim);
  for l=1:30
      par0=((S'*W0*S)^(-1))*S'*W0*temp;
      w0=sum(W0);
      var0=(w0*((temp-S*par0).^2))/(sum(w0));
      for h=1:indim
          W0(h,h)=exp(-0.5*taw*((temp(h)-S(h,:)*par0)^2)/var0);
      end
  end 
  
  
      I11=zeros(ttt);
      for h=1:indim
          I11=I11+((W0(h,h)^2)*((temp(h)-S(h,:)*par0)^2)*H(h,:)'*H(h,:))/indim;
      end
      
      
      I22=zeros(ppp);
      for h=1:indim
          I22=I22+((W0(h,h)^2)*((temp(h)-S(h,:)*par0)^2)*S(h,:)'*S(h,:))/indim;
      end
  
      I12=zeros(ttt,ppp);
      for h=1:indim
      I12=I12+((W0(h,h)^2)*((temp(h)-S(h,:)*par0)^2)*H(h,:)'*S(h,:))/indim;
      end
        
     der=zeros(ttt,1);
     for h=1:indim
          der=der+W0(h,h)*(temp(h)-S(h,:)*par0)*H(h,:)';
     end
      
T1=der'*(I11-I12*(I22^(-1))*I12')^(-1)*der;

TraoAlpha(i)=T1/indim;
 

 %%%%%%%  Huber Rao%%%%
tun=1.34;
  W0=eye(indim);
  for l=1:30
      par0=((S'*W0*S)^(-1))*S'*W0*temp;
      w0=sum(W0);
      var0=(w0*((temp-S*par0).^2))/(sum(w0));
      for h=1:indim
          dh=((temp(h)-S(h,:)*par0))/sqrt(var0);
          W0(h,h)=min([1 tun/abs(dh)]);
      end
  end 
  
  
      I11=zeros(ttt);
      for h=1:indim
          I11=I11+((W0(h,h)^2)*((temp(h)-S(h,:)*par0)^2)*H(h,:)'*H(h,:))/indim;
      end
      
      
      I22=zeros(ppp);
      for h=1:indim
          I22=I22+((W0(h,h)^2)*((temp(h)-S(h,:)*par0)^2)*S(h,:)'*S(h,:))/indim;
      end
  
      I12=zeros(ttt,ppp);
      for h=1:indim
      I12=I12+((W0(h,h)^2)*((temp(h)-S(h,:)*par0)^2)*H(h,:)'*S(h,:))/indim;
      end
        
     der=zeros(ttt,1);
     for h=1:indim
          der=der+W0(h,h)*(temp(h)-S(h,:)*par0)*H(h,:)';
     end
      
T1=der'*(I11-I12*(I22^(-1))*I12')^(-1)*der;

THuber(i)=T1/indim;

  
%%%%%%% Hypothesis 1 Wald Robust %%%%%%%%%%%%%%%%%   
  W1=eye(indim);
  for l=1:30
      par1=((Htilda'*W1*Htilda)^(-1))*Htilda'*W1*temp;
      w1=sum(W1);
      var1=(w1*((temp-Htilda*par1).^2))/(sum(w1));
      for h=1:indim
          W1(h,h)=exp(-0.5*taw*((temp(h)-Htilda(h,:)*par1)^2)/var1);
      end
  end 
W2=eye(indim);
      for h=1:indim
          W2(h,h)=1-taw*((temp(h)-Htilda(h,:)*par1)^2)/var1;
      end

L=(1/(var1*indim))*(Htilda'*(W1.^2)*Htilda);
K=(1/(var1*indim))*(Htilda'*(W1.*W2)*Htilda);
V=(K^(-1))*L*(K^(-1));
whitenedpar=(V^(-0.5))*(par1);
tarpar=whitenedpar(1:ttt);
TwaldAlpha(i)=indim*(tarpar)'*(tarpar);
end

end