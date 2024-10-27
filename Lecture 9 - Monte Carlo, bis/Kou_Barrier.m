% Price a Down&Out Barrier Call Option
% Kou model
%% Parameters
K=100; S0=105; D=90; T=1; r=0.01;
M=round(52*T); % weekly monitoring
Nsim=1e6;
sigma=0.4; p=0.7; lambdap=7; lambdam=8; lambda=3;
char_exp=@(u) -sigma^2/2*u.^2+1i*u*lambda.*...
    (p./(lambdap-1i*u)-(1-p)./(lambdam+1i*u));
drift=r-char_exp(-1i); dt=T/M;
%% Standard MC
X=zeros(Nsim,M+1);
Z=randn(Nsim,M); NJ=poissrnd(lambda*dt,Nsim,M);
for i=1:M
    X(:,i+1)=X(:,i)+drift*dt+sigma*sqrt(dt)*Z(:,i);
    for j=1:Nsim
        for ii=1:NJ(j,i)
            if rand<p
                X(j,i+1)=X(j,i+1)+exprnd(1/lambdap);
            else
                X(j,i+1)=X(j,i+1)-exprnd(1/lambdam);
            end
        end
    end
end
S=S0*exp(X); 
disp('r-n check: must be zero')
[check,~,CI_check]=normfit(S(:,end)-S0*exp(r*T))
[price,~,CI]=normfit( exp(-r*T)*max(S(:,end)-K,0).*...
    (min(S,[],2)>D) )
%% AV MC (only Brownian)
X=zeros(Nsim,M+1); XAV=X;
Z=randn(Nsim,M); NJ=poissrnd(lambda*dt,Nsim,M);
for i=1:M
    X(:,i+1)=X(:,i)+drift*dt+sigma*sqrt(dt)*Z(:,i);
    XAV(:,i+1)=XAV(:,i)+drift*dt-sigma*sqrt(dt)*Z(:,i);
    for j=1:Nsim
        for ii=1:NJ(j,i)
            if rand<p
                JJ=exprnd(1/lambdap);
            else
                JJ=-exprnd(1/lambdam);
            end
            X(j,i+1)=X(j,i+1)+JJ;
            XAV(j,i+1)=XAV(j,i+1)+JJ;            
        end
    end
end
S=S0*exp(X); SAV=S0*exp(XAV);
[price,~,CI]=normfit( exp(-r*T)*(...
    max(S(:,end)-K,0).*(min(S,[],2)>D) +...
    max(SAV(:,end)-K,0).*(min(SAV,[],2)>D) )/2 )
%% AV MC (also jump intensity)
X=zeros(Nsim,M+1); XAV=X;
Z=randn(Nsim,M); 
U=rand(Nsim,M);
NJ=icdf('Poisson',U,lambda*dt);
NJAV=icdf('Poisson',1-U,lambda*dt);
for i=1:M
    X(:,i+1)=X(:,i)+drift*dt+sigma*sqrt(dt)*Z(:,i);
    XAV(:,i+1)=XAV(:,i)+drift*dt-sigma*sqrt(dt)*Z(:,i);
    for j=1:Nsim
        for ii=1:NJ(j,i)
            if rand<p
                JJ=exprnd(1/lambdap);
            else
                JJ=-exprnd(1/lambdam);
            end
            X(j,i+1)=X(j,i+1)+JJ;
        end
        for ii=1:NJAV(j,i)
            if rand<p
                JJ=exprnd(1/lambdap);
            else
                JJ=-exprnd(1/lambdam);
            end
            XAV(j,i+1)=XAV(j,i+1)+JJ;            
        end
    end
end
S=S0*exp(X); SAV=S0*exp(XAV);
[price,~,CI]=normfit( exp(-r*T)*(...
    max(S(:,end)-K,0).*(min(S,[],2)>D) +...
    max(SAV(:,end)-K,0).*(min(SAV,[],2)>D) )/2 )

%% AV MC (also jump size)
X=zeros(Nsim,M+1); XAV=X;
Z=randn(Nsim,M); NJ=poissrnd(lambda*dt,Nsim,M);
for i=1:M
    X(:,i+1)=X(:,i)+drift*dt+sigma*sqrt(dt)*Z(:,i);
    XAV(:,i+1)=XAV(:,i)+drift*dt-sigma*sqrt(dt)*Z(:,i);
    for j=1:Nsim
        for ii=1:NJ(j,i)
            if rand<p
                u=rand;
                JJ=icdf('Exponential',u,(1/lambdap));
                JJAV=icdf('Exponential',1-u,(1/lambdap));
            else
                u=rand;
                JJ=-icdf('Exponential',u,(1/lambdam));
                JJAV=-icdf('Exponential',1-u,(1/lambdam));
            end
            X(j,i+1)=X(j,i+1)+JJ;
            XAV(j,i+1)=XAV(j,i+1)+JJAV;            
        end
    end
end
S=S0*exp(X); SAV=S0*exp(XAV);
[price,~,CI]=normfit( exp(-r*T)*(...
    max(S(:,end)-K,0).*(min(S,[],2)>D) +...
    max(SAV(:,end)-K,0).*(min(SAV,[],2)>D) )/2 )
%% CV MC
Ef=S0*exp(r*T);
% Sample f=S(T) -> estimate alpha
Nsim2=1000;
X=zeros(Nsim2,M+1);
Z=randn(Nsim2,M); NJ=poissrnd(lambda*dt,Nsim2,M);
for i=1:M
    X(:,i+1)=X(:,i)+drift*dt+sigma*sqrt(dt)*Z(:,i);
    for j=1:Nsim2
        for ii=1:NJ(j,i)
            if rand<p
                X(j,i+1)=X(j,i+1)+exprnd(1/lambdap);
            else
                X(j,i+1)=X(j,i+1)-exprnd(1/lambdam);
            end
        end
    end
end
S=S0*exp(X); 
f=S(:,end); g=exp(-r*T)*max( f-K,0 ).*(min(S,[],2)>D);
VC=cov(f,g); alpha=-VC(1,2)/VC(1,1);
% Compute the price
X=zeros(Nsim,M+1);
Z=randn(Nsim,M); NJ=poissrnd(lambda*dt,Nsim,M);
for i=1:M
    X(:,i+1)=X(:,i)+drift*dt+sigma*sqrt(dt)*Z(:,i);
    for j=1:Nsim
        for ii=1:NJ(j,i)
            if rand<p
                X(j,i+1)=X(j,i+1)+exprnd(1/lambdap);
            else
                X(j,i+1)=X(j,i+1)-exprnd(1/lambdam);
            end
        end
    end
end
S=S0*exp(X); 
f=S(:,end); g=exp(-r*T)*max( f-K,0 ).*(min(S,[],2)>D);
[price,~,CI]=normfit( g+alpha*(f-Ef) )
