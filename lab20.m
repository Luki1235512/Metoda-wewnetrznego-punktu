clc
close all
clear

A = [3 2 1; 5 4 3; 1 1 0];
b = [120; 300; 50];
c = [12; 8; 10];
c = -c;

wynik = fkcja(A,b,c,10^(-8))

function x = fkcja(A,b,c,tol)
m = size(A, 1);
A = [A eye(m)];
c(end+1:end+m)=0;
n = size(A, 2);

delta = 0.02;
sigma = 0.9;

x = ones(n,1);
y = ones(m,1);
z = ones(n,1);
    
while 1
    rp = A*x-b;
    rd = A'*y +z-c;
    xz = dot(x,z);
    
    if  norm(rp) < tol && norm(rd) < tol && xz < tol
        break 
    end
    
    my = (1/n)*xz;
    tau = my*delta;
    X = diag(x);
    Z = diag(z);
    e = ones(length(X),1);

    delta_y = (A*(X*(Z^(-1))*A'))\(-rp+A*(x-tau*(Z^(-1))*e-X*((Z^(-1))*rd)));
    delta_z  = -A'*delta_y-rd;
    delta_x = (-X*Z^(-1)*delta_z+tau*Z^(-1)*e-x);
    
    alpha_p = 1;
    alpha_d = 1;
    for j = 1:n
        if x(j) + alpha_p * delta_x(j) < 0
            alpha_p = -x(j)/delta_x(j);  
        end
        if z(j) + alpha_d * delta_z(j) < 0
            alpha_d = -z(j)/delta_z(j);
        end
    end
    
    alpha_p = alpha_p * sigma;
    alpha_d = alpha_d * sigma;
    
    x = x + alpha_p * delta_x;
    y = y + alpha_d * delta_y;
    z = z + alpha_d * delta_z;
    
end
x = x(1:end-m);    
end


