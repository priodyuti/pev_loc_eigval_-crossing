%%
clear all
clc

path = 'wheel_RR_deloc_graph_1.txt'

%path = 'ER_599901.txt'

fd = fopen(path,'rt');

formatSpec = '%d %d';
sizeA = [2 Inf];

Y = fscanf(fd, formatSpec, sizeA);
Y = Y';
fclose(fd);

N = max(max(Y(:,1)),max(Y(:,2)));
A = zeros(N, N);
t = size(Y, 1);

for i=1:t
    A(Y(i,1),Y(i,2)) = 1;
    A(Y(i,2),Y(i,1)) = 1;
end    
clear Y;

deg = sum(A,2);
[Max_Deg Max_Index] = max(deg);
total_deg = sum(A(:));

fprintf('Number of Nodes: %d\n', N);
fprintf('Number of edges: %d\n',total_deg/2 ); 
fprintf('Max_Deg: %d\n', Max_Deg); 
fprintf('Max_Index: %d\n', Max_Index); 
fprintf('Average_deg: %d\n', total_deg/N); 

x = ones(N,1);

mu = 0.5;
f = 2.6;
I = eye(N,N);
L = 18;

M = f*(1-mu)*I + (f*mu/(3*L))*A;

t = 300000;
k = 1500;
m = ceil(t/k);
data = zeros(m+1, N);
d = norm(x,2);
x = x./d;
data(1,:) = x';
j = 2;
for i = 1:t
   x = M*x;
   d = norm(x,2);
   x = x./d;
   if mod(i,k) == 1
    %fprintf('i = %d\n',i);
    data(j,:) = x';
       %i = i + k;
    j = j + 1;
   end          
end
%display(x)
lambda = x'*M*x;
%save(['stable_state_loc.mat'],'x');
save('wheel_random_500_deloc.mat','data');

IPR = sum(x.^4);
fprintf('IPR = %0.16f max_eigval = %f\n', IPR,lambda);

[evec,eval1] = eig(M);
eigval = diag(eval1);
Ev = evec(:,N);
%d = norm(Ev,1);  %normalized with norm 1
%Ev = Ev./d;
IPR = sum(Ev.^4);
fprintf('IPR = %0.16f max_eigval = %f sec_max_eigval = %f\n', IPR,eigval(N),eigval(N-1));

clear all


