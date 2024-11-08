close all
clear
% Knock&Out Call on a basket of two underlying assets (not correlated)
T=1;
S10=2; S20=1; K=0.9;
D1=0.8; U1=4; D2=0.5; U2=2;
r=0.1/100; sigma1=0.6; sigma2=0.4;
x1min=log(D1/S10); x1max=log(U1/S10);
x2min=log(D2/S20); x2max=log(U2/S20);
Mt=100; dt=T/Mt;
N1=50; N2=50;
x1=linspace(x1min,x1max,N1+1); dx1=x1(2)-x1(1);
x2=linspace(x2min,x2max,N2+1); dx2=x2(2)-x2(1);
num_points=(N1+1)*(N2+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               i+1
%   i-(N2-1)     i    i+(N2+1)
%               i-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundaries
W=[1:N2+1];
S=[1:N2+1:num_points];
N=[(N2+1):N2+1:num_points];
E=[S(end)+(0:N2)];
Boundary=sort(unique([N E S W]));
% Matrix
M=spalloc(num_points,num_points,num_points*5);
for i=1:num_points
    if min( abs(i-Boundary) )==0 % i is a boundary index
        M(i,i)=1;
    else
        % time derivative
        M(i,i)=-1/dt;
        % no derivative term
        M(i,i)=M(i,i)-r;
        % first order derivative w.r.t x1
        M(i,i+(N2+1))= (r-sigma1^2/2)/(2*dx1);
        M(i,i-(N2+1))= -(r-sigma1^2/2)/(2*dx1);
        % first order derivative w.r.t x2
        M(i,i+1)= (r-sigma2^2/2)/(2*dx2);
        M(i,i-1)= -(r-sigma2^2/2)/(2*dx2);
        % second order derivative w.r.t x1
        M(i,i+(N2+1))=M(i,i+(N2+1))+(sigma1^2/2)/dx1^2;
        M(i,i)=M(i,i)-2*(sigma1^2/2)/dx1^2;
        M(i,i-(N2+1))=M(i,i-(N2+1))+(sigma1^2/2)/dx1^2;
        % second order derivative w.r.t x2
        M(i,i+1)=M(i,i+1)+(sigma2^2/2)/dx2^2;
        M(i,i)=M(i,i)-2*(sigma2^2/2)/dx2^2;
        M(i,i-1)=M(i,i-1)+(sigma2^2/2)/dx2^2;
    end
end
% Backward in time
[X1,X2]=meshgrid(x1,x2);
V=max( 0.5*(S10*exp(X1(:))+S20*exp(X2(:))) - K,0);
V(Boundary)=0;
for j=Mt:-1:1
    rhs=-V/dt; 
    rhs(Boundary)=0; % Boundary Condition
    V=M\rhs;
end
Vmat=reshape(V,size(X1));
figure
surf(X1,X2,Vmat)
price=griddata(S10*exp(X1), S20*exp(X2), Vmat, S10, S20 )












