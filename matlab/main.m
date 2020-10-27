% ***************************************************************************
% Polynomial Interpolation
% Code Library
% ***************************************************************************
% Author: Chaobin
% Email:  chaubyZou@163.com
% Date: October 2020
% ***************************************************************************
% Language: Matlab
% Also available in: Python
% Required library: None
% ***************************************************************************

clear;
close all;

%% Specify the data point
% q_given = [0, 1.6, 3.2, 2, 4, 0.2, 1.2; % position
%             0, 0, 0, 0, 0, 0, 0; % zero velocity
%            0, 0, 0, 0, 0, 0, 0]'; % zero acceleration

q_given = [0, 1.6, 3.2, 2, 4, 0.2, 1.2; % position
           0, 1.0, 2.0, -2.0, -1.0, 0, 0; % varying velocity
           0, 1, 2, 3, 2, 1, 0]'; % varying acceleration

t_given = [0, 1, 3, 4.5, 6, 8, 10]'; % time

% time for interpolation
t = t_given(1):0.01:t_given(end);


%% ************************ Linear interpolation *******************************
linear_interpolation = LinearInterpolation('Linear', q_given, t_given);
linear_trajectory = zeros(length(t), 3); % N x 3 array: position, velocity, acceleration

for i = 1:length(t)
    linear_trajectory(i,:) = linear_interpolation.getPosition(t(i));
end

figure('position', [200, 200, 800, 600])
subplot(3,1,1)
hold on
plot(t_given, q_given(:, 1), 'ro')
plot(t, linear_trajectory(:,1), 'k')
hold off
grid on
title('Linear interpolation')
xlabel('time')
ylabel('position')
xlim([t_given(1)-1, t_given(end)+1])
ylim([min(q_given(:,1)) - 1, max(q_given(:,1)) + 1])

subplot(3,1,2)
hold on
plot(t_given, q_given(:, 2), 'ro')
plot(t, linear_trajectory(:,2), 'k')
hold off
grid on
xlabel('time')
ylabel('velocity')
xlim([t_given(1)-1, t_given(end)+1])

subplot(3,1,3)
hold on
plot(t_given, q_given(:, 3), 'ro')
plot(t, linear_trajectory(:,3), 'k')
hold off
grid on
xlabel('time')
ylabel('acceleration')
xlim([t_given(1)-1, t_given(end)+1])


%% *********************** Parabolic interpolation *****************************
parabolic_interpolation = ParabolicInterpolation('Parabolic', q_given, t_given);
parabolic_trajectory = zeros(length(t), 3); % N x 3 array: position, velocity, acceleration

for i = 1:length(t)
    parabolic_trajectory(i,:) = parabolic_interpolation.getPosition(t(i));
end

figure('position', [230, 200, 800, 600])
subplot(3,1,1)
hold on
plot(t_given, q_given(:, 1), 'ro')
plot(t, parabolic_trajectory(:,1), 'k')
hold off
grid on
title('Parabolic interpolation')
xlabel('time')
ylabel('position')
xlim([t_given(1)-1, t_given(end)+1])
ylim([min(q_given(:,1)) - 1, max(q_given(:,1)) + 1])

subplot(3,1,2)
hold on
plot(t_given, q_given(:, 2), 'ro')
plot(t, parabolic_trajectory(:,2), 'k')
hold off
grid on
xlabel('time')
ylabel('velocity')
xlim([t_given(1)-1, t_given(end)+1])

subplot(3,1,3)
hold on
plot(t_given, q_given(:, 3), 'ro')
plot(t, parabolic_trajectory(:,3), 'k')
hold off
grid on
xlabel('time')
ylabel('acceleration')
xlim([t_given(1)-1, t_given(end)+1])


%% ************************ Cubic interpolation ******************************
cubic_interpolation = CubicInterpolation('Cubic', q_given, t_given);
cubic_trajectory = zeros(length(t), 3); % N x 3 array: position, velocity, acceleration

for i = 1:length(t)
    cubic_trajectory(i,:) = cubic_interpolation.getPosition(t(i));
end

figure('position', [260, 200, 800, 600])
subplot(3,1,1)
hold on
plot(t_given, q_given(:, 1), 'ro')
plot(t, cubic_trajectory(:,1), 'k')
hold off
grid on
title('Cubic interpolation')
xlabel('time')
ylabel('position')
xlim([t_given(1)-1, t_given(end)+1])
ylim([min(q_given(:,1)) - 1, max(q_given(:,1)) + 1])

subplot(3,1,2)
hold on
plot(t_given, q_given(:, 2), 'ro')
plot(t, cubic_trajectory(:,2), 'k')
hold off
grid on
xlabel('time')
ylabel('velocity')
xlim([t_given(1)-1, t_given(end)+1])

subplot(3,1,3)
hold on
plot(t_given, q_given(:, 3), 'ro')
plot(t, cubic_trajectory(:,3), 'k')
hold off
grid on
xlabel('time')
ylabel('acceleration')
xlim([t_given(1)-1, t_given(end)+1])


%% *********************** Polynomial 5 interpolation ****************************
polynomial5_interpolation = Polynomial5Interpolation('Polynomial 5', q_given, t_given);
polynomial5_trajectory = zeros(length(t), 3); % N x 3 array: position, velocity, acceleration

for i = 1:length(t)
    polynomial5_trajectory(i,:) = polynomial5_interpolation.getPosition(t(i));
end

figure('position', [290, 200, 800, 600])
subplot(3,1,1)
hold on
plot(t_given, q_given(:, 1), 'ro')
plot(t, polynomial5_trajectory(:,1), 'k')
hold off
grid on
title('Polynomial of degree 5 interpolation')
xlabel('time')
ylabel('position')
xlim([t_given(1)-1, t_given(end)+1])
ylim([min(q_given(:,1)) - 1, max(q_given(:,1)) + 1])

subplot(3,1,2)
hold on
plot(t_given, q_given(:, 2), 'ro')
plot(t, polynomial5_trajectory(:,2), 'k')
hold off
grid on
xlabel('time')
ylabel('velocity')
xlim([t_given(1)-1, t_given(end)+1])

subplot(3,1,3)
hold on
plot(t_given, q_given(:, 3), 'ro')
plot(t, polynomial5_trajectory(:,3), 'k')
hold off
grid on
xlabel('time')
ylabel('acceleration')
xlim([t_given(1)-1, t_given(end)+1])


%% ***************************** Comparison ***********************************
figure('position', [350, 200, 800, 600])
subplot(3,1,1)
hold on
plot(t_given, q_given(:, 1), 'ro')
plot(t, linear_trajectory(:,1), 'k')
plot(t, parabolic_trajectory(:,1), 'b')
plot(t, cubic_trajectory(:,1), 'g')
plot(t, polynomial5_trajectory(:,1), 'm')
hold off
grid on
title('Comparison')
xlabel('time')
ylabel('position')
legend('given pos', 'linear', 'parabolic', 'cubic', 'poly 5')
xlim([t_given(1)-1, t_given(end)+3])
ylim([min(q_given(:,1)) - 1, max(q_given(:,1)) + 1])

subplot(3,1,2)
hold on
plot(t_given, q_given(:,2), 'ro')
plot(t, linear_trajectory(:,2), 'k')
plot(t, parabolic_trajectory(:,2), 'b')
plot(t, cubic_trajectory(:,2), 'g')
plot(t, polynomial5_trajectory(:,2), 'm')
hold off
grid on
xlabel('time')
ylabel('velocity')
legend('given vel', 'linear', 'parabolic', 'cubic', 'poly 5')
xlim([t_given(1)-1, t_given(end)+3])

subplot(3,1,3)
hold on
plot(t_given, q_given(:,3), 'ro')
plot(t, linear_trajectory(:,3), 'k')
plot(t, parabolic_trajectory(:,3), 'b')
plot(t, cubic_trajectory(:,3), 'g')
plot(t, polynomial5_trajectory(:,3), 'm')
hold off
grid on
xlabel('time')
ylabel('acceleration')
legend('given acc', 'linear', 'parabolic', 'cubic', 'poly 5')
xlim([t_given(1)-1, t_given(end)+3])
