%% =========== beta = 0.5 omega = 0.6 ========================
% ============= 平面波相 =============================
% close all
DENA1 = dlmread('DENA1.DAT');
DENA2 = dlmread('DENA2.DAT');
DENB1 = dlmread('DENB1.DAT');
DENB2 = dlmread('DENB2.DAT');
PHASEA1 = dlmread('PHASEA1.DAT');
PHASEA2 = dlmread('PHASEA2.DAT');
PHASEB1 = dlmread('PHASEB1.DAT');
PHASEB2 = dlmread('PHASEB2.DAT');
x = linspace(-100, 100, 513);

figure(1);
% multiply with 100
plot(x, 100*DENA1, 'r-', x, 100*DENA2, 'b:', 'LineWidth', 2);
xlabel('x/k_0^{-1}');
ylabel('\rho(10^{-2})');
legend('\rho_\uparrow', '\rho_\downarrow');

%%
% close all
figure(2);
plot(x, PHASEA1, 'r-', x, PHASEA2, 'b:', 'LineWidth', 0.8);
xlabel('x/k_0^{-1}');
% ylabel('\phi');
legend('\phi_\uparrow', '\phi_\downarrow');

%% 
figure(3);
% multiply with 100
plot(x, 100*DENB1, 'r-', x, 100*DENB2, 'b:', 'LineWidth', 2);
xlabel('x/k_0^{-1}');
ylabel('\rho(10^{-2})');
legend('\rho_\uparrow', '\rho_\downarrow');

%% 
figure(4);
plot(x, PHASEB1, 'r-', x, PHASEB2, 'b:', 'LineWidth', 0.8);
xlabel('x/k_0^{-1}');
% ylabel('\phi');
legend('\phi_\uparrow', '\phi_\downarrow');