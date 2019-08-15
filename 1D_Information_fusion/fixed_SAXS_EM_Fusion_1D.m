function [C_theta, C_theta_noa] = fixed_SAXS_EM_Fusion_1D(num, angle, c_theta)
% SAXS-EM fusion in 1D
% Notations are changing in the draft, please refer to the following
% C = \hat\rho 
% C_hat = \tilde\rho
% M = F

fprintf('This is a 2D SAXS-EM Fusion Checking Program.\n');

%% Input from user
% Define theta 2n*1 matrix,C(theta) 2n*1 matrix
theta_half = (angle.*(pi/180));
theta = [theta_half ; theta_half+pi];
C = [c_theta;c_theta];

% Define avaerage density at the circle with radius = 1:
fun = @(x) ((2.*(sin(x)).^5 + cos(x) + 2.*(cos(x)).^3).^2).^2;
I_radius1 = integral(fun,0,2*pi);

%% Variable declaration
%Set up matrix V(column matrix with 1):
V = ones(num*2,1);

%Set matrix M and M_negative1(take out the zeroth)
colmultiplier = [-num:-1,1:num];
M = exp(1i*theta*colmultiplier);
M_negative1=M^(-1);
% Set up for unknown a
% syms a real
syms a

% Execution
eqn = M_negative1 * (C - a*V);
result_a = solve(a^2 + eqn' * eqn - I_radius1/(2*pi),a);
result_a_simplified = result_a;

%% Get C_hat(a)(which is \tilde\rho' in draft) (equation 10)
result_a_simplified_negative1 = result_a_simplified';
for iteration1 = 1:length(result_a_simplified_negative1)
    C_hat_a{iteration1} = M_negative1 * (C - result_a_simplified_negative1(iteration1)*V);
end

%% Get C_(theta) (based on a) (equation 4 for any theta value, interpolated)
syms the
numMultiplier = (-num):(num);
for ite=1:length(result_a_simplified_negative1)
    C_hat{ite} = [C_hat_a{ite}(1:num,1);result_a_simplified(ite);C_hat_a{ite}(num+1:2*num,1)];
    C_theta{ite} = sum( C_hat{ite}.*exp(1i*(numMultiplier)'*the));
    C_theta{ite} = vpa(C_theta{ite},4);
end

%% Solve C_theta_noa without a

M_1 = exp(1i*theta*numMultiplier); %2n*(2n+1) matrix

M_1(2*num+1,:) = zeros(1,2*num+1);
M_1(2*num+1,1) = -1;
M_1(2*num+1,2*num+1) = 1;
C_noa = [C;0]; % by setting the last item in C is 0

C_hat_k_solved = M_1^(-1) * C_noa;

%% Get C_theta without a (equation 4)

C_hat_k_cell = C_hat_k_solved; 
C_theta_noa = sum( C_hat_k_cell.*exp(1i*(numMultiplier)'*the));

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
legend('y1','y2','y3','y4');
end
