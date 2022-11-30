% ============= 色散关系 =============
close all;
clear all;
kx = linspace(-2, 2, 100);
kr = 1;
m = 1;
E_recoil = kr^2/(2*m);
delta = 0.5;
% omega = 2*E_recoil;

for i = 1:3
    figure;
%     omega = i*E_recoil;
    omega = 2*i*E_recoil;
    % energy
    E_up = (kx.^2 + kr^2)/(2*m) + sqrt(kx.^2*kr^2/m^2 + kx.*kr*delta/m + delta^2/4 + omega^2/4);
    E_down = (kx.^2 + kr^2)/(2*m) - sqrt(kx.^2*kr^2/m^2 + kx.*kr*delta/m + delta^2/4 + omega^2/4);
    plot(kx, E_up, '--', kx, E_down, 'LineWidth', 1);
    xlabel('{\it k_x}/{\it k_r}','FontName','Times New Roman');
    ylabel('{\it E}/{\it E_r}', 'FontName','Times New Roman');
    axis([-2 2 -2, 5]);
    disp(i);
end
