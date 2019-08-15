function [C_theta, C_theta_noa] = SAXS_EM_Fusion_1D(num, angle, c_theta)
% SAXS-EM fusion in 1D
% Notations are changing in the draft, please refer to the following
% C = \hat\rho 
% C_hat = \tilde\rho
% M = F

fprintf('This is a 2D SAXS-EM Fusion Checking Program.\n');

%% Input from user
% Define theta
theta = zeros(1,num);
for w = 1:num
    theta(w) = angle(w) * (pi/180);
end
e_min = num + 1;
e_max = num * 2;
for e = e_min:e_max
    theta(e) = theta(e-num) + pi;
end
theta = theta'; %theta is a 2n*1 matrix

% Define c(theta)(also known as C):
for e = e_min:e_max
    c_theta(e) = c_theta(e-num);
end
c_theta = c_theta'; %c_theta is a 2n*1 matrix
C = c_theta;

% Define avaerage density at the circle with radius = 1:
fun = @(x) ((2.*(sin(x)).^5 + cos(x) + 2.*(cos(x)).^3).^2).^2;
I_radius1 = integral(fun,0,2*pi);

%% Variable declaration
%Set up matrix V(column matrix with 1):
V = ones(num*2,1);

%Set matrix M and M_negative1(take out the zeroth)
M = zeros(num*2);
colmultiplier = (-num):(num);
colmultiplier((length(colmultiplier) + 1)/2) = [];
for m = 1:(num*2)
    for n = 1:(num*2)
        M(m,n) = exp(1i*colmultiplier(n)*theta(m));
    end
end
M_negative1 = M^(-1);

% Set up for unknown a
% syms a real
syms a

% Execution
eqn = M_negative1 * (C - a*V);
result_a = solve(a^2 + eqn' * eqn - I_radius1/(2*pi),a);
result_a_simplified = result_a;

%% Get C_hat(a)(which is \tilde\rho' in draft) (equation 10)
result_a_simplified_negative1 = result_a_simplified';
iteration1 = 1;
while iteration1 <= length(result_a_simplified_negative1)
    C_hat_a{iteration1} = M_negative1 * (C - result_a_simplified_negative1(iteration1)*V);
    iteration1 = iteration1 + 1;
end

%% Get C_(theta) (based on a) (equation 8 for any theta value, interpolated)
syms the
ite = 1;
for i = 1:length(result_a_simplified_negative1)
    C_theta{i} = 0;
end
while ite <= length(result_a_simplified_negative1)
    C_hat{ite} = [C_hat_a{ite}(1:num,1);result_a_simplified(ite);C_hat_a{ite}(num+1:2*num,1)];
    j = -1*num;
    for i = 1:2*num+1
        C_theta{ite} = C_theta{ite} + C_hat{ite}(i)*exp(1i*j*the);
        j = j+1;
    end
    C_theta{ite} = vpa(C_theta{ite},4);
    ite = ite + 1;
end

%% Solve C_theta_noa without a
numMultiplier = (-num):(num);
for j = 1:(num*2)
    for k = 1:(num*2 + 1)
        M_1(j,k) = exp(1i*numMultiplier(k)*theta(j)); %2n*(2n+1) matrix
    end
end 

M_1(2*num+1,:) = zeros(1,2*num+1);
M_1(2*num+1,1) = -1;
M_1(2*num+1,2*num+1) = 1;
C_noa = [C;0]; % by setting the last item in C is 0

X = sym ('x', [2*num + 1 1]);
C_hat_k_solved = M_1^(-1) * C_noa;

%% Get C_theta without a (equation 4)

C_hat_k_cell = C_hat_k_solved; 
C_theta_noa = 0;
j = -1*num;
for i = 1:(num*2+1)
    C_theta_noa = C_theta_noa + C_hat_k_cell(i,1)*exp(1i*j*the);
    j = j+1;
end

%% Print out
% print out conditions:
fprintf('I(1) is: %f.\n',I_radius1);
for s = 1:length(result_a_simplified)
    fprintf('%f.\n',result_a_simplified(s));
end

figure
% plot
x = 0:pi/20:2*pi;
y1 = subs(C_theta{1,1}, x);
y2 = subs(C_theta{1,2}, x);
y3 = subs(C_theta_noa, x);
y4 = (2.*(sin(x)).^5 + cos(x) + 2.*(cos(x)).^3).^2;
plot(x,y1,'-.r','LineWidth',3); hold on;
plot(x,y2,'--g','LineWidth',3);
plot(x,y3,':k','LineWidth',1.5);hold on;
plot(x,y4,'--bs','LineWidth',1.5,'MarkerEdgeColor','k','MarkerSize',3);
set(gcf,'color','w');

end
