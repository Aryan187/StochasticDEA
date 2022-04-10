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
tau = zeros(N,M+H+S,M+H+S);
tau(1,:,:) = table2array(readtable('tau1.csv'));
tau(2,:,:) = table2array(readtable('tau2.csv'));
tau(3,:,:) = table2array(readtable('tau3.csv'));
tau(4,:,:) = table2array(readtable('tau4.csv'));
tau(5,:,:) = table2array(readtable('tau5.csv'));
tau(6,:,:) = table2array(readtable('tau6.csv'));
tau(7,:,:) = table2array(readtable('tau7.csv'));
tau(8,:,:) = table2array(readtable('tau8.csv'));
tau(9,:,:) = table2array(readtable('tau9.csv'));
tau(10,:,:) = table2array(readtable('tau10.csv'));
tau(11,:,:) = table2array(readtable('tau11.csv'));
tau(12,:,:) = table2array(readtable('tau12.csv'));
tau(13,:,:) = table2array(readtable('tau13.csv'));
tau(14,:,:) = table2array(readtable('tau14.csv'));
tau(15,:,:) = table2array(readtable('tau15.csv'));
tau(16,:,:) = table2array(readtable('tau16.csv'));
tau(17,:,:) = table2array(readtable('tau17.csv'));
tau(18,:,:) = table2array(readtable('tau18.csv'));
tau(19,:,:) = table2array(readtable('tau19.csv'));
tau(20,:,:) = table2array(readtable('tau20.csv'));
map = [1,2,3,4,7,8,9,10,5,6,11,12];
tmp = tau;
for i = 1:N
    for j = 1:12
        for k = 1:12
            tau(i,map(j),map(k)) = tmp(i,j,k);
        end
    end
end
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
alpha = 0.025;
negpinv = -1.96;
effs = [];
%Stochastic DEA with copula
for k = 1:N
    %The variables are X = (l1,l2,...ln,u1,u2,....,un,theta)
    %min C.T*X
    c = zeros(2*N+1,1);
    c(2*N+1) = 1;
    %||A_(i)*x-b_(i)|| <= d(i)*x - g(i) for each constraint
    id = 1;
    for i = 1:M1
        for l = i:M1
            A_i = zeros(2*N+1,2*N+1);
            b_i = zeros(2*N+1,1);
            d_i = zeros(2*N+1,1);
            g_i = 0;
            for j = 1:N
                A_i(j,j) = tau(j,i,l)*s_x(l,j);
                d_i(j) = tau(j,i,l)*u_x(l,j);
            end
            A_i(k,2*N+1) = -1*tau(k,i,l)*s_x(l,k);
            d_i(2*N+1) = -1*tau(k,i,l)*u_x(l,k);
            d_i = d_i/negpinv;
            socConstraints(id) = secondordercone(A_i,b_i,d_i,g_i);
            id = id + 1;
        end
    end
    %||A_(i)*x-b_(i)|| <= d(i)*x - g(i) for each constraint
    for i = M1+1:M
        for l = i:M
            A_i = zeros(2*N+1,2*N+1);
            b_i = zeros(2*N+1,1);
            d_i = zeros(2*N+1,1);
            g_i = 0;
            for j = 1:N
                A_i(N+j,N+j) = tau(j,i,l)*s_x(l,j);
                d_i(N+j) = tau(j,i,l)*u_x(l,j);
            end
            A_i(N+k,2*N+1) = -1*tau(k,i,l)*s_x(l,k);
            d_i(2*N+1) = -1*tau(k,i,l)*u_x(l,k);
            d_i = d_i/negpinv;
            socConstraints(id) = secondordercone(A_i,b_i,d_i,g_i);
            id = id + 1;
        end
    end
    %||A_(i)*x-b_(i)|| <= d(i)*x - g(i) for each constraint
    for i = 1:H
        for l = i:i
            A_i = zeros(2*N+1,2*N+1);
            b_i = zeros(2*N+1,1);
            d_i = zeros(2*N+1,1);
            g_i = 0;
            for j = 1:N
                A_i(j,j) = -1*tau(j,M+i,M+l)*s_z(l,j);
                A_i(j,N+j) = tau(j,M+i,M+l)*s_z(l,j);
                d_i(j) = -1*tau(j,M+i,M+l)*u_z(l,j);
                d_i(N+j) = tau(j,M+i,M+l)*u_z(l,j);
            end
            d_i = d_i/negpinv;
            socConstraints(id) = secondordercone(A_i,b_i,d_i,g_i);
            id = id + 1;
        end
    end
    %||A_(i)*x-b_(i)|| <= d(i)*x - g(i) for each constraint
    for i = 1:S1
        A_i = zeros(2*N+1,2*N+1);
        b_i = zeros(2*N+1,1);
        d_i = zeros(2*N+1,1);
        g_i = 0;
        for j = 1:N
            A_i(j,j) = tau(j,M+H+i,M+H+l)*s_y(i,j);
            d_i(j) = -1*tau(j,M+H+i,M+H+l)*u_y(i,j);
        end
        b_i(k) = tau(k,M+H+i,M+H+l)*s_y(i,k);
        g_i = -1*tau(k,M+H+i,M+H+l)*u_y(i,k)/negpinv;
        d_i = d_i/negpinv;
        socConstraints(id) = secondordercone(A_i,b_i,d_i,g_i);
        id = id + 1;
    end
    %||A_(i)*x-b_(i)|| <= d(i)*x - g(i) for each constraint
    for i = S1+1:S
        for l = i:S
            A_i = zeros(2*N+1,2*N+1);
            b_i = zeros(2*N+1,1);
            d_i = zeros(2*N+1,1);
            g_i = 0;
            for j = 1:N
                A_i(N+j,N+j) = tau(j,M+H+i,M+H+l)*s_y(l,j);
                d_i(N+j) = -1*tau(j,M+H+i,M+H+l)*u_y(l,j);
            end
            b_i(N+k) = tau(k,M+H+i,M+H+l)*s_y(l,k);
            g_i = -1*tau(k,M+H+i,M+H+l)*u_y(l,k)/negpinv;
            d_i = d_i/negpinv;
            socConstraints(id) = secondordercone(A_i,b_i,d_i,g_i);
            id = id + 1;
        end
    end
    lb = zeros(2*N+1,1);
    ub = Inf(2*N+1,1);
    %display(id);
    x = coneprog(c,socConstraints,[],[],[],[],lb,ub,options);
    effs = [effs x(2*N+1)];
end

display(effs);
