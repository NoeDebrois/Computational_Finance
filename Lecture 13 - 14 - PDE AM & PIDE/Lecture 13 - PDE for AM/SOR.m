function xnew=SOR(A,rhs,x)
xnew=zeros(size(x)); 
n=length(A);
maxiter=100;
tol=1e-7;
w=1.5;
for i=1:maxiter
    for j=1:n
        y=(rhs(j)-A(j,1:j-1)*xnew(1:j-1,1)...
            -A(j,j+1:end)*x(j+1:end,1))/A(j,j);
        xnew(j)=x(j)+w*(y-x(j));
    end
    if norm(xnew-x)<tol
        break
    else
        x=xnew;
    end
end