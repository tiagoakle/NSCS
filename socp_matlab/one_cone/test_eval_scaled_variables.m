clear all
%Test the eval scaled variables method 

%Choose a problem size 
n = 10;
x = zeros(n,1);
s = zeros(n,1);
max_norm = 10; %Select a bound for x0

%Choose a random x,s
%Pick a pair of random points inside the unit sphere of 
%size one

d = randn(n-1,1); %Pick a direction
d = d/norm(d);
x(1)   = max_norm*rand(1);
x(2:n) = x(1)*d*sqrt(rand(1)); % Pick a distance from the origin

d = randn(n-1,1); %Pick a direction
d = d/norm(d);
s(1)   = max_norm*rand(1);
s(2:n) = s(1)*d*sqrt(rand(1)); % Pick a distance from the origin

%Call the method 
[H,W,Wi,lambda] = eval_scaled_variables(n,x,s);

%Tests
fprintf('|Hx-s|           %g \n',norm(H*x-s));
fprintf('|W*W-H|          %g \n',norm(W*W-H));
fprintf('|W*x-lambda|     %g \n',norm(W*x-lambda));
fprintf('|Wi*s-lambda|    %g \n',norm(Wi*s-lambda));
fprintf('|Wi*H*Wi-eye(n)| %g \n',norm(Wi*H*Wi-eye(n)));

