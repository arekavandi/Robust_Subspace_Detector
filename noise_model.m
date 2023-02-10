function out=noise_model(N,dim,ep,sigma1,mu1,sigma2,mu2,sw)
n=N*dim;

pn=floor((1-ep)*n);
sn=n-pn;

noise1=random('normal',mu1,sigma1,pn,1);
noise2=random('normal',mu2,sigma2,sn,1);
if sw==1;
    noise2=abs(noise2);
end

noise=[noise1;noise2];
out=shuffle(noise);


out=reshape(out,dim,N);


 function v=shuffle(v)
     v=v(randperm(length(v)));
 end
end