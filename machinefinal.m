%% 1. Parameters
VL = 6600; S = 2500000; f = 50; p = 1;
Ra = 0.07; Xl = 0.9;
Pmech = 20000; Pcore = 12000;
slots = 60; cs = 4; Nf = 100;

Mdata = [5000 10000 15000 20000 25000 30000];   % AT per pole
Edata = [1500 2800 3700 4200 4400 4600];        % V per phase

%% 2. Basic calculations
Vph = VL / sqrt(3);
Ia = S / (sqrt(3) * VL);
ns = 120 * f / (2 * p);
ws = 2 * pi * ns / 60;

q = slots / (2 * p * 3);
alpha = 60 / q;
kd = sind(q * alpha / 2) / (q * sind(alpha / 2));
kw = kd;

Ma_FL = 1.35 * cs * q * Ia * kw;

%% 3. OCC and air gap lines
M = 5000:30000;
E = spline(Mdata, Edata, M);
M = [0 M];
E = [0 E];
If = 2 * M / Nf;

m1 = Edata(1) / Mdata(1);
M0 = interp1(E, M, Vph);
m2 = Vph / M0;
If0 = 2 * M0 / Nf;

%% 4. Short circuit characteristic
Ersc = Ia * Xl;
Mrsc = Ersc / m1;
Mfsc = Mrsc + Ma_FL;
Ifsc = 2 * Mfsc / Nf;
Isc = Ia * M / Mfsc;
Isc0 = Ia * If0 / Ifsc;

%% 5. Plot
figure('Color', 'white');

yyaxis left
plot(If, E, 'b-', 'LineWidth', 2);
hold on;
plot(If, m1 * M, 'k--', 'LineWidth', 1.5);
plot(If, m2 * M, 'm-.', 'LineWidth', 1.5);
ylabel('Open Circuit Phase Voltage (V)');
ylim([0 10000]);

yyaxis right
plot(If, Isc, 'r-', 'LineWidth', 2);
ylabel('Short Circuit Current Isc (A)');
ylim([0 1000]);

xlabel('Field Current If (A)');
title('Synchronous Generator Characteristics (OCC, Air-Gap Line, SCC)');
legend('OCC', 'Air-Gap Line', 'Modified Air-Gap Line', 'SCC', 'Location', 'northwest');
grid on;

%% 6. Synchronous reactance
Xu1 = m1 * M0 / Isc0;
Xs1 = Vph / Isc0;

Ifu = 2 * (Vph / m1) / Nf;
SCRu = Ifu / Ifsc;
SCRs = If0 / Ifsc;
Zbase = Vph / Ia;
Xu2 = Zbase / SCRu;
Xs2 = Zbase / SCRs;

fprintf('--- Synchronous Reactance ---\n');
fprintf('Method 1 (OCC/SCC): Xu = %.2f ohm | Xs = %.2f ohm\n', Xu1, Xs1);
fprintf('Method 2 (SCR)    : Xu = %.2f ohm | Xs = %.2f ohm\n', Xu2, Xs2);

%% 7. Full load and half load performance
pfs = [0.8 0.9 1.0 0.9 0.8];
names = {'0.8 Lag', '0.9 Lag', 'Unity', '0.9 Lead', '0.8 Lead'};
loads = [1.0 0.5];

for n = 1:length(loads)
    L = loads(n);
    I = Ia * L;

    fprintf('\n====== Load = %.0f%% ======\n', L * 100);

    for m = 1:length(pfs)
        pf = pfs(m);

        if m <= 3
            phi = acosd(pf);
        else
            phi = -acosd(pf);
        end

        Ic = I * (cosd(phi) - 1i * sind(phi));

        % exact (ampere-turn) method
        Er = Vph + Ic * (Ra + 1i * Xl);
        Mr = interp1(E, M, abs(Er));
        Ma = Ma_FL * L;
        beta = phi + angle(Er) * 180 / pi;
        Mf = sqrt(Mr^2 + Ma^2 + 2 * Mr * Ma * sind(beta));
        If_ex = 2 * Mf / Nf;
        E0_ex = sqrt(3) * interp1(M, E, Mf);
        VR_ex = 100 * (E0_ex - VL) / VL;

        % approximate method
        Ef = abs(Vph + Ic * (Ra + 1i * Xs1));
        If_ap = 2 * Ef / (Nf * m2);
        E0_ap = sqrt(3) * Ef;
        VR_ap = 100 * (E0_ap - VL) / VL;

        Sout = S * L;
        P = Sout * pf;
        Q = Sout * sind(phi);

        Pcu = 3 * I^2 * Ra;
        Pin = P + Pcu + Pcore + Pmech;
        T = Pin / ws;

        fprintf('Case: %s | S = %.2f kVA\n', names{m}, Sout / 1000);
        fprintf('  Exact : If = %.1f A | E0 = %.1f V | VR = %.2f %%\n', If_ex, E0_ex, VR_ex);
        fprintf('  Approx: If = %.1f A | E0 = %.1f V | VR = %.2f %%\n', If_ap, E0_ap, VR_ap);
        fprintf('  P = %.2f kW | Q = %.2f kVAR\n', P / 1000, Q / 1000);
        fprintf('  Pin = %.2f kW | T = %.2f N.m\n', Pin / 1000, T);
        fprintf('----------------------------\n');
    end
end