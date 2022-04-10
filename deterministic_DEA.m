clear all;
M = 6;
M1 = 4;
H = 3;
S = 3;
S1 = 1;
N = 20;
data = table2array(readtable('data.xlsx'));
u_x = rand(M,N)*100;
u_y = rand(S,N)*100;
u_z = rand(H,N)*100;
s_x = rand(M,N);
s_y = rand(S,N);
s_z = rand(H,N);
alpha = 0.025;
negpinv = -1.96;
effs = [];
for i = 1:N
    u_x(1,i) = data(i,2);
    u_x(2,i) = data(i,4);
    u_x(3,i) = data(i,6);
    u_x(4,i) = data(i,8);
    u_z(1,i) = data(i,10);
    u_z(2,i) = data(i,12);
    u_z(3,i) = data(i,14);
    u_y(1,i) = data(i,16);
    u_x(5,i) = data(i,18);
    u_x(6,i) = data(i,20);
    u_y(2,i) = data(i,22);
    u_y(3,i) = data(i,24);
end
options = optimoptions('linprog','Display','none');
%Deterministic DEA
for k = 1:N
    %The variables are X = (l1,l2,...ln,u1,u2,....,un,theta)
    %min C.T*X
    c = zeros(2*N+1,1);
    c(2*N+1) = 1;
    %AX <= b
    b = zeros(M+H+S,1);
    A = zeros(M+H+S,2*N+1);
    for i = 1:M1
        A(i,1:N) = u_x(i,1:N);
        A(i,2*N+1) = -1*u_x(i,k);
    end
    for i = M1+1:M
        A(i,N+1:2*N) = u_x(i,1:N);
        A(i,2*N+1) = -1*u_x(i,k);
    end
    for i = 1:H
        A(M+i,1:N) = -1*u_z(i,1:N);
        A(M+i,N+1:2*N) = u_z(i,1:N);
    end
    for i = 1:S1
        A(M+H+i,1:N) = -1*u_y(i,1:N);
        b(M+H+i) = -1*u_y(i,k);
    end
    for i = S1+1:S
        A(M+H+i,N+1:2*N) = -1*u_y(i,1:N);
        b(M+H+i) = -1*u_y(i,k);
    end
    lb = zeros(2*N+1,1);
    ub = Inf(2*N+1,1);
    x = linprog(c,A,b,[],[],lb,ub,options);
    effs = [effs x(2*N+1)];
end
display(effs);