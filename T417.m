clc;
clear all;
close all;
rho = 1;
c = 2.5e6;
k_prime = 2;
Lh = 100e6;
Tm = 0;
k = k_prime / (rho * c);
L = 1;
days = 50;
total_time = 3600 * 24 * days;
dx = 0.125;
dt = 3600;
x = 0 : dx : L;
nx = length(x);
nt = round(total_time / dt);
T_initial = 2.0;
T_boundary = -10;
if k * (dt / (dx^2)) < 0.5
    fprintf("Stable solution\n");
else
    fprintf("Unstable solution\n");
end
T = T_initial * ones(1, nx);
T(1) = T_boundary;
H = c * T;
for i = 1 : nx
    if T(i) >= Tm
    H(i) = H(i) + Lh;
    end
end
i_p = [];
i_t = [];
T_h = {};
t_h = [];
H_t = c * Tm + 0.5 * Lh;
alpha = k_prime * dt / (rho * dx^2);
for n = 1:nt
    c_t = (n-1) * dt;
    H_old = H;
    T_old = zeros(1, nx);
    for i = 1 : nx
        if H_old(i) < c * Tm
            T_old(i) = H_old(i) / c;
        elseif H_old(i) > c * Tm + Lh
            T_old(i) = (H_old(i) - Lh) / c;
        else
            T_old(i) = Tm;
        end
    end
    H_new = H_old;
    for i = 2 : nx-1
        H_new(i) = H_old(i) + alpha * (T_old(i+1) - 2*T_old(i) + T_old(i-1));
    end
    H_new(1) = c * T_boundary;
    H_new(end) = H_new(end-1);
    for i = 2 : nx-1
        if (H_old(i) > H_t) && (H_new(i) <= H_t)
            xcap = (H_t - H_old(i)) / (H_new(i) - H_old(i));
            e_t = c_t + xcap * dt;
            i_p(end+1) = x(i);
            i_t(end+1) = e_t;
        end
    end
    H = H_new;
    for i = 1 : nx
        if H(i) < c * Tm
        T(i) = H(i) / c;
        elseif H(i) > c * Tm + Lh
        T(i) = (H(i) - Lh) / c;
        else
        T(i) = Tm;
        end
    end
    if (mod(n, 200) == 0) || (n == nt)
        T_h{end+1} = T;
        t_h(end+1) = c_t + dt;
    end
end
lambda = 0.3224;
x_in_analy = 2 * lambda * sqrt(k * total_time);
T_analytical = zeros(1, nx);
for i = 1 : nx
    if x(i) <= x_in_analy
        T_analytical(i) = T_boundary + (Tm - T_boundary) * (erf(x(i) / (2 * sqrt(k * total_time))) / erf(lambda));
    else
        T_analytical(i) = T_initial - (T_initial - Tm) * (erfc(x(i) / (2 * sqrt(k * total_time))) / erfc(lambda));
    end
end
target = 0.5;
min_diff = L;
target_idx = 1;
for i = 1 : length(x)
    current_diff = abs(x(i) - target);
    if current_diff < min_diff
        min_diff = current_diff;
        target_idx = i;
    end
end
target_x = x(target_idx);
figure();
plot(i_t/3600/24, i_p * 100, 'bo', 'MarkerSize', 5, 'LineWidth', 1.5);
hold on;
t_ana = linspace(0, total_time, 100);
X_ana = 2 * lambda * sqrt(k .* t_ana);
plot(t_ana/3600/24, X_ana * 100, 'r-', 'LineWidth', 2);
xlabel('Time, days');
ylabel('Position of phase change boundary, cm');
title('Movement of the Phase Change Boundary');
legend('Numerical solution', 'Analytic', 'Location', 'SouthEast');
grid on;
figure();
monitor_temp_history = zeros(1, length(t_h));
for i = 1 : length(t_h)
    T_snap = T_h{i};
    monitor_temp_history(i) = T_snap(target_idx);
end
plot(t_h/3600/24, monitor_temp_history, 'bo', 'MarkerSize', 5, 'LineWidth', 1.5);
hold on;
T_history_ana = ones(size(t_ana));
X_history_ana = 2 * lambda * sqrt(k .* t_ana);
for i = 1 : length(t_ana)
    if t_ana(i) == 0
        T_history_ana(i) = T_initial;
    else
        if target_x > X_history_ana(i)
        T_history_ana(i) = T_initial - (T_initial - Tm) * (erfc(target_x / (2 * sqrt(k * t_ana(i)))) / erfc(lambda));
        else
        T_history_ana(i) = T_boundary + (Tm - T_boundary) * (erf(target_x / (2 * sqrt(k * t_ana(i)))) / erf(lambda));
        end
    end
end
plot(t_ana/3600/24, T_history_ana, 'r-', 'LineWidth', 2);
xlabel('Time, days');
ylabel('Temperature, °C');
title(sprintf('Temperature History at x=%.0f cm', target_x * 100));
legend('Numerical solution', 'Analytic', 'Location', 'best');
ylim([T_boundary, T_initial]);
grid on;
