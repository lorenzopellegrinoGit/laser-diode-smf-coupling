lambda_laser = 0.8;      % Lunghezza d'onda in um
n_clad = 3.1710;      
n_co=3.3067;
             
k0 = 2*pi / lambda_laser;
%verifica che siano effettivamente single mode
d=0.4; %thickness of the core in order to have single mode dopo calcoli su carta con solita formula

beta_vector = linspace(k0*n_clad, k0*n_co, 5000);
kappa = sqrt(k0^2 * n_co^2 - beta_vector.^2);
gamma = sqrt(beta_vector.^2 - k0^2 * n_clad^2);

first_member = tan(kappa * d);
second_member = (2 * kappa .* gamma) ./ (kappa.^2 - gamma.^2);


figure('Name', 'Ricerca Beta_0 Guida SM');
plot(beta_vector, first_member, 'r', 'LineWidth', 1.5); hold on;
plot(beta_vector, second_member, 'b', 'LineWidth', 1.5);
grid on;
ylim([-15, 15]);
xlabel('\beta (\mum^{-1})');
ylabel('Valore delle funzioni');
legend('First member: \tan(\kappa d)', 'Second member: 2\kappa\gamma / (\kappa^2 - \gamma^2)');

beta0=25.5729; %um-1

kappa0 = sqrt(k0^2 * n_co^2 - beta0^2);
gamma0 = sqrt(beta0^2 - k0^2 * n_clad^2);

y= linspace(-3, 3, 5000); %visto che il core va da 0 a -d 
A = 1; 
% B/A = -delta / kappa. 
% la guida è simmetrica quindi l'esponenziale n3 è uguale a quello nel substrate
delta0 = gamma0;
B = -A * (delta0 / kappa0);
E0 = zeros(size(y));

top_idx = y >= 0;
core_idx  = y < 0 & y >= -d;
sub_idx   = y < -d;

E0(top_idx) = A * exp(-delta0 * y(top_idx));
E0(core_idx)  = A * cos(kappa0 * y(core_idx)) + B * sin(kappa0 * y(core_idx));
E0(sub_idx)   = (A * cos(kappa0 * d) - B * sin(kappa0 * d)) * exp(gamma0 * (y(sub_idx) + d));


figure('Name', 'Profilo del Campo');
plot(y, E0, 'g' , 'LineWidth', 2);
hold on;
xline(0, 'k--'); xline(-d, 'k--');
grid on;
xlabel('Asse trasversale y (\mum)');
ylabel('Ampiezza del Campo');

dy = y(2) - y(1);

% %integrale di |E0|^2 calcolato come somma delle altezze per la base
Z0 = 376.73; % Impedenza del vuoto in Ohm
neff_laser = beta0 / k0; % Indice efficace del modo

% Integrale della potenza con ammettenza modale Y = neff / Z0
P0 = (neff_laser / Z0) * sum(abs(E0).^2) * dy;
E0_norm = E0 / sqrt(P0); % normalizzazione


figure('Name', 'Profilo del Campo Normalizzato (SM)');
plot(y, E0_norm,'m', 'LineWidth', 2);
hold on;
xline(0, 'k--'); xline(-d, 'k--'); % Confini del core
grid on;
title('Profilo del Campo Normalizzato E_{0,norm}(y)');
xlabel('Asse trasversale y [\mum]');
ylabel('Ampiezza del Campo normalizzata');

% PARAMETRI DI BASE
w0_x = 25;      % um (Waist orizzontale, 50 um / 2)
MFD = 5;        % um (Filtro fibra)
w_f = MFD / 2;  % um (Waist della fibra)
z_target_1s = [100, 1000]; % Distanze in um (100 um e 1 mm)

% 1: FITTING DEL MODO VERTICALE (Per trovare w0_y)
% Il modo Gaussian ideale è E_gauss = exp(-(y/w0_y)^2)
% Trovo w0_y misurando dove il campo E0 scende a 1/e del suo massimo
% (oppure l'intensità scende a 1/e^2)

[max_E, max_idx] = max(E0_norm);
y_center = y(max_idx);
E_target = max_E * exp(-1); % Valore 1/e
% Trova l'indice dove il campo attraversa il valore 1/e
idx_1e = find(E0_norm(max_idx:end) <= E_target, 1) + max_idx - 1;
w0_y = abs(y(idx_1e) - y_center);

fprintf('Waist verticale fittato: w0_y = %.4f um\n', w0_y);


% 2: PROPAGAZIONE & OVERLAP ANALITICO

disp('--- RISULTATI DI ACCOPPIAMENTO');

for i = 1:length(z_target_1s)  %La funzione length() calcola quanti elementi ci sono nell'array z_target_1s cioe 2
    % La variabile i è l'indice (o contatore) che, ad ogni giro del ciclo, si incrementa di 1
    % Estrae l'$i$-esimo valore dall'array e lo assegna alla variabile temporanea $z$.
    z = z_target_1s(i);
    
    % Asse X (Orizzontale) formule nel report
    zR_x = (pi * w0_x^2) / lambda_laser;
    wx_z = w0_x * sqrt(1 + (z/zR_x)^2);
    Rx_z = z * (1 + (zR_x/z)^2);
    
    % Asse Y (Verticale) formule nel report
    zR_y = (pi * w0_y^2) / lambda_laser;
    wy_z = w0_y * sqrt(1 + (z/zR_y)^2);
    Ry_z = z * (1 + (zR_y/z)^2);
    
    % Calcolo Efficienza Asse X
    term1_x = (wx_z / w_f) + (w_f / wx_z);
    term2_x = (pi * wx_z * w_f) / (lambda_laser * Rx_z);
    eta_x = 4 / (term1_x^2 + term2_x^2);
    
    % Calcolo Efficienza Asse Y
    term1_y = (wy_z / w_f) + (w_f / wy_z);
    term2_y = (pi * wy_z * w_f) / (lambda_laser * Ry_z);
    eta_y = 4 / (term1_y^2 + term2_y^2);
    
    % Efficienza Totale
    eta_tot = eta_x * eta_y;
    Coupling_Loss_dB = -10 * log10(eta_tot);
    
    fprintf('Distanza z = %d um:\n', z);
    fprintf('  w_x(z) = %.2f um, R_x(z) = %.2f um -> eta_x = %.2f%%\n', wx_z, Rx_z, eta_x*100);
    fprintf('  w_y(z) = %.2f um, R_y(z) = %.2f um -> eta_y = %.2f%%\n', wy_z, Ry_z, eta_y*100);
    fprintf('  Efficienza Totale = %.4f%% | Coupling Loss = %.2f dB\n\n', eta_tot*100, Coupling_Loss_dB);
end

% Grafico 1: Coupling Loss vs Distanza 
z_array = linspace(1, 1000, 500); % da 1 a 1000 um
eta_tot_array = zeros(size(z_array));

for i = 1:length(z_array)
    z_curr = z_array(i);
    
    % Asse X
    wx = w0_x * sqrt(1 + (z_curr/((pi * w0_x^2)/lambda_laser))^2);
    Rx = z_curr * (1 + (((pi * w0_x^2)/lambda_laser)/z_curr)^2);
    term1_x = (wx / w_f) + (w_f / wx);
    term2_x = (pi * wx * w_f) / (lambda_laser * Rx);
    eta_x_curr = 4 / (term1_x^2 + term2_x^2);
    
    % Asse Y
    wy = w0_y * sqrt(1 + (z_curr/((pi * w0_y^2)/lambda_laser))^2);
    Ry = z_curr * (1 + (((pi * w0_y^2)/lambda_laser)/z_curr)^2);
    term1_y = (wy / w_f) + (w_f / wy);
    term2_y = (pi * wy * w_f) / (lambda_laser * Ry);
    eta_y_curr = 4 / (term1_y^2 + term2_y^2);
    
    % Efficienza Totale (pura sovrapposizione geometrica)
    eta_tot_array(i) = eta_x_curr * eta_y_curr;
end

Loss_dB_array = -10 * log10(eta_tot_array);

figure('Name', 'Coupling Loss vs Distance');
plot(z_array, Loss_dB_array, 'r', 'LineWidth', 2);
grid on;
xlabel('Air Gap Distance z [\mum]', 'FontWeight', 'bold');
ylabel('Coupling Loss [dB]', 'FontWeight', 'bold');
title('Coupling Loss Evolution in Free-Space (Analytical Model)', 'FontWeight', 'bold');
% Evidenzia i due punti calcolati nel loop precedente
hold on;
idx_100 = find(z_array>=100, 1);
plot(100, Loss_dB_array(idx_100), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
plot(1000, Loss_dB_array(end), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
legend('Coupling Loss', 'Case 1: z = 100 \mum', 'Case 2: z = 1 mm', 'Location', 'southeast');


theta = linspace(0, 2*pi, 100);
x_fiber = w_f * cos(theta); 
y_fiber = w_f * sin(theta);

figure('Name', 'Laser beam vs fiber core', 'Position', [100 100 1000 400]);

% Subplot 1: z = 100 um
subplot(1,2,1);
wx_100 = w0_x * sqrt(1 + (100/zR_x)^2);
wy_100 = w0_y * sqrt(1 + (100/zR_y)^2);
plot(wx_100 * cos(theta), wy_100 * sin(theta), 'r--', 'LineWidth', 1.5); hold on;
plot(x_fiber, y_fiber, 'b-', 'LineWidth', 1.5);
fill(x_fiber, y_fiber, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
grid on; axis equal;
title(' z = 100 \mum', 'FontWeight', 'bold');
xlabel('X [\mum]'); ylabel('Y [\mum]');
legend('Laser Beam', 'Fiber Acceptance');
xlim([-150 150]); ylim([-150 150]); 

% Subplot 2: z = 1000 um
subplot(1,2,2);
wx_1000 = w0_x * sqrt(1 + (1000/zR_x)^2);
wy_1000 = w0_y * sqrt(1 + (1000/zR_y)^2);
plot(wx_1000 * cos(theta), wy_1000 * sin(theta), 'g--', 'LineWidth', 1.5); hold on;
plot(x_fiber, y_fiber, 'b-', 'LineWidth', 1.5);
fill(x_fiber, y_fiber, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
grid on; axis equal;
title('z = 1000 \mum', 'FontWeight', 'bold');
xlabel('X [\mum]'); ylabel('Y [\mum]');
legend('Laser Beam', 'Fiber Acceptance');
xlim([-1000 1000]); ylim([-1000 1000]);

% 4: ABCD MATRICES, CLIPPING & AUTOMATION

disp('--- RISULTATI ACCOPPIAMENTO CON LENTI');

% a. Scelta lenti commerciali Thorlabs
f1 = 15000;   % LA1509-B
D1 = 12700;   % Diametro 1/2 pollice in um
f2 = 150000;  % LA1417-B
D2 = 25400;   % Diametro 1 pollice in um

d1 = f1;      % Distanza Laser-L1
d2 = 50000;   % Wc (Volo libero)
d3 = f2;      % Distanza L2-Fibra

% b. Distanze di Rayleigh iniziali
zR_x = (pi * w0_x^2) / lambda_laser;
zR_y = (pi * w0_y^2) / lambda_laser;

% c. Calcolo dell'ingombro su L1 e perdita per Clipping
w1_x = w0_x * sqrt(1 + (d1/zR_x)^2);
w1_y = w0_y * sqrt(1 + (d1/zR_y)^2);

% La trasmissione per l'apertura circolare di L1 (domina l'asse Y)
a1 = D1 / 2; % Raggio della lente
T_clip_L1 = 1 - exp(-2 * (a1 / w1_y)^2);
Loss_clip_dB = -10 * log10(T_clip_L1);

% d. Matrici ABCD (Dall'ingresso a L2 fino alla fibra)
M_d1 = [1, d1; 0, 1];
M_L1 = [1, 0; -1/f1, 1];
M_d2 = [1, d2; 0, 1];
M_L2 = [1, 0; -1/f2, 1];
M_d3 = [1, d3; 0, 1];

M_tot = M_d3 * M_L2 * M_d2 * M_L1 * M_d1;
A = M_tot(1,1); B = M_tot(1,2); C = M_tot(2,1); D = M_tot(2,2);

% Parametri complessi q iniziali
q0_x = 1j * zR_x;
q0_y = 1j * zR_y;

% Trasformazione ABCD
q3_x = (A * q0_x + B) / (C * q0_x + D);
q3_y = (A * q0_y + B) / (C * q0_y + D);

inv_q3_x = 1 / q3_x;
inv_q3_y = 1 / q3_y;

w3_x = sqrt(-lambda_laser / (pi * imag(inv_q3_x))); %im(invq3) è zr nuova dopo la matrice
R3_x = 1 / real(inv_q3_x);
w3_y = sqrt(-lambda_laser / (pi * imag(inv_q3_y)));
R3_y = 1 / real(inv_q3_y);

% e. Efficienze Finali
term1_x = (w3_x / w_f) + (w_f / w3_x);
term2_x = (pi * w3_x * w_f) / (lambda_laser * R3_x);
eta_x_lens = 4 / (term1_x^2 + term2_x^2);

term1_y = (w3_y / w_f) + (w_f / w3_y);
term2_y = (pi * w3_y * w_f) / (lambda_laser * R3_y);
eta_y_lens = 4 / (term1_y^2 + term2_y^2);

eta_tot_geom = eta_x_lens * eta_y_lens;
Loss_geom_dB = -10 * log10(eta_tot_geom);

Total_System_Loss = Loss_geom_dB + Loss_clip_dB;

fprintf('Ingombro su L1: w1_x = %.1f um, w1_y = %.1f um\n', w1_x, w1_y);
fprintf('Trasmissione L1 (Clipping) = %.1f%% (%.2f dB di perdita)\n', T_clip_L1*100, Loss_clip_dB);
fprintf('Parametri alla fibra: w3_x = %.2f um, w3_y = %.2f um\n', w3_x, w3_y);
fprintf('Parametri alla fibra: R3_x = %.2f, R3_y = %.2f um\n', R3_x, R3_y);
fprintf('Efficienza geometrica (X) = %.2f%%\n', eta_x_lens * 100);
fprintf('Efficienza geometrica (Y) = %.2f%%\n', eta_y_lens * 100);
fprintf('Coupling Loss Geometrica = %.2f dB\n', Loss_geom_dB);
fprintf('Coupling Loss TOTALE (incluso clipping) = %.2f dB\n\n', Total_System_Loss);

% --- Grafico 3: Impronta Spaziale alla Faccetta della Fibra (DOPO LE LENTI)
theta = linspace(0, 2*pi, 200);

% Ellisse del Laser dopo le lenti (Usa i parametri calcolati con matrici ABCD)
x_laser_lens = w3_x * cos(theta);
y_laser_lens = w3_y * sin(theta);

% Cerchio della Fibra
x_fiber = w_f * cos(theta);
y_fiber = w_f * sin(theta);

figure('Name', 'Laser beam vs fiber core after lenses');
% Disegna il fascio laser trasformato
plot(x_laser_lens, y_laser_lens, 'r--', 'LineWidth', 2); hold on;

% Disegna e colora la fibra
plot(x_fiber, y_fiber, 'b-', 'LineWidth', 2);
fill(x_fiber, y_fiber, 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none'); 

grid on; axis equal;
xlabel('Horizontal Axis X [\mum]', 'FontWeight', 'bold');
ylabel('Vertical Axis Y [\mum]', 'FontWeight', 'bold');
title('Laser beam vs fiber core after lenses', 'FontWeight', 'bold');

xlim([-250 250]);
ylim([-10 10]); 

