clear all;
M = 6;
M1 = 4;
H = 3;
S = 3;
S1 = 1;
N = 20;
data = table2array(readtable('data.xlsx'));
%display(data);
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
    s_x(1,i) = data(i,3);
    u_x(2,i) = data(i,4);
    s_x(2,i) = data(i,5);
    u_x(3,i) = data(i,6);
    s_x(3,i) = data(i,7);
    u_x(4,i) = data(i,8);
    s_x(4,i) = data(i,9);
    u_z(1,i) = data(i,10);
    s_z(1,i) = data(i,11);
    u_z(2,i) = data(i,12);
    s_z(2,i) = data(i,13);
    u_z(3,i) = data(i,14);
    s_z(3,i) = data(i,15);
    u_y(1,i) = data(i,16);
    s_y(1,i) = data(i,17);
    u_x(5,i) = data(i,18);
    s_x(5,i) = data(i,19);
    u_x(6,i) = data(i,20);
    s_x(6,i) = data(i,21);
    u_y(2,i) = data(i,22);
    s_y(2,i) = data(i,23);
    u_y(3,i) = data(i,24);
    s_y(3,i) = data(i,25);
end
options = optimoptions('coneprog','Display','none');
%Stochastic DEA
for k = 1:N
    %The variables are X = (l1,l2,...ln,u1,u2,....,un,theta)
    %min C.T*X
    c = zeros(2*N+1,1);
    c(2*N+1) = 1;
    %||A_(i)*x-b_(i)|| <= d(i)*x - g(i) for each constraint
    for i = 1:M1
        A_i = zeros(2*N+1,2*N+1);
        b_i = zeros(2*N+1,1);
        d_i = zeros(2*N+1,1);
        g_i = 0;
        for j = 1:N
            A_i(j,j) = s_x(i,j);
            d_i(j) = u_x(i,j);
        end
        A_i(k,2*N+1) = -1*s_x(i,k);
        d_i(2*N+1) = -1*u_x(i,k);
        d_i = d_i/negpinv;
        socConstraints(i) = secondordercone(A_i,b_i,d_i,g_i);
    end
    %||A_(i)*x-b_(i)|| <= d(i)*x - g(i) for each constraint
    for i = M1+1:M
        A_i = zeros(2*N+1,2*N+1);
        b_i = zeros(2*N+1,1);
        d_i = zeros(2*N+1,1);
        g_i = 0;
        for j = 1:N
            A_i(N+j,N+j) = s_x(i,j);
            d_i(N+j) = u_x(i,j);
        end
        A_i(N+k,2*N+1) = -1*s_x(i,k);
        d_i(2*N+1) = -1*u_x(i,k);
        d_i = d_i/negpinv;
        socConstraints(i) = secondordercone(A_i,b_i,d_i,g_i);
    end
    %||A_(i)*x-b_(i)|| <= d(i)*x - g(i) for each constraint
    for i = 1:H
        A_i = zeros(2*N+1,2*N+1);
        b_i = zeros(2*N+1,1);
        d_i = zeros(2*N+1,1);
        g_i = 0;
        for j = 1:N
            A_i(j,j) = -1*s_z(i,j);
            A_i(j,N+j) = s_z(i,j);
            d_i(j) = -1*u_z(i,j);
            d_i(N+j) = u_z(i,j);
        end
        d_i = d_i/negpinv;
        socConstraints(M+i) = secondordercone(A_i,b_i,d_i,g_i);
    end
    %||A_(i)*x-b_(i)|| <= d(i)*x - g(i) for each constraint
    for i = 1:S1
        A_i = zeros(2*N+1,2*N+1);
        b_i = zeros(2*N+1,1);
        d_i = zeros(2*N+1,1);
        g_i = 0;
        for j = 1:N
            A_i(j,j) = s_y(i,j);
            d_i(j) = -1*u_y(i,j);
        end
        b_i(k) = s_y(i,k);
        g_i = -1*u_y(i,k)/negpinv;
        d_i = d_i/negpinv;
        socConstraints(M+H+i) = secondordercone(A_i,b_i,d_i,g_i);
    end
    %||A_(i)*x-b_(i)|| <= d(i)*x - g(i) for each constraint
    for i = S1+1:S
        A_i = zeros(2*N+1,2*N+1);
        b_i = zeros(2*N+1,1);
        d_i = zeros(2*N+1,1);
        g_i = 0;
        for j = 1:N
            A_i(N+j,N+j) = s_y(i,j);
            d_i(N+j) = -1*u_y(i,j);
        end
        b_i(N+k) = s_y(i,k);
        g_i = -1*u_y(i,k)/negpinv;
        d_i = d_i/negpinv;
        socConstraints(M+H+i) = secondordercone(A_i,b_i,d_i,g_i);
    end
    lb = zeros(2*N+1,1);
    ub = Inf(2*N+1,1);
    x = coneprog(c,socConstraints,[],[],[],[],lb,ub);
    effs = [effs x(2*N+1)];
end

display(effs);
