function dydt = Kinetic_Equation(~,y,Para_set)
% ODEs of the model

% Variables

G = y(1);                                                                  % Glucose
X = y(2);                                                                  % Xylose
E = y(3);                                                                  % Ethanol
Nx = y(4);  Ng = y(5);                                                     % Cell densities (OD)
Rx = y(6);  Rg = y(7);                                                     % Precursor
O = y(8);                                                                  % Oxygen
Cx = y(9);                                                                 % NAD+ in YX

% Parameters
a_x = Para_set(1);   K_x = Para_set(2);   k_x = Para_set(5);               % Xylose consumption of YX strain
a_g = Para_set(3);   K_g = Para_set(4);                                    % Glucose consumption of YG strain

a_epx = Para_set(6); K_epx = Para_set(7);                                  % Ethanol production of YX
a_epg = Para_set(6); K_epg = Para_set(7);                                  % Ethanol production of YG
a_ecx = Para_set(8); K_ecx = Para_set(9); k_ecx = Para_set(10);            % Ethanol reassimilation of YX
a_ecg = Para_set(8); K_ecg = Para_set(9); k_ecg = Para_set(11);            % Ethanol reassimilation of YG
a_nx = Para_set(12); K_nx = Para_set(13);                                  % Flux to Cell growth in YX
a_ng = Para_set(12); K_ng = Para_set(13);                                  % Flux to Cell growth in YG

a_mx = Para_set(14);  a_mg = Para_set(14);  K_mx = 1e-6;  K_mg = 1e-6;       % Flux to Maintenance 


b_x = Para_set(15);   b_g = Para_set(16);                                  % yielding constant from sugar to precursor
b_epx = Para_set(17); b_epg = Para_set(17);                                % yielding constant from precursor to ethanol
b_ecx = Para_set(18); b_ecg = Para_set(18);                                % yielding constant from ethanol to precursor
g_x = Para_set(19);  g_g = Para_set(19);                                   % conversion coefficient from precursor to biomass

v_b = Para_set(20);                                                        % oxygen supplement
a_o = Para_set(21);   K_o = Para_set(22);                                  % oxygen consumption rate
b_o = Para_set(23);                                                        % convertion coefficient from oxygen to NAD+
k_xc = Para_set(24);                                                       % NAD+ consumed in utilizing unit weight xylose
K_C = Para_set(25);                                                        % half-velocity NAD concentration for xylose consumption

a_bc = Para_set(26);  K_bc = Para_set(27);                                 % Basal NAD+ consumption

% Flux rates
J_g   = Hill(a_g,  K_g,  G);                                               % glucose consumption rate
J_x   = Hill(a_x,  K_x,  X) /(1+k_x*E) * Cx/(Cx+K_C);                      % xylose consumption rate
J_ng  = Hill(a_ng, K_ng, Rg) ;                                             % Flux rate to growth (YG)
J_nx  = Hill(a_nx, K_nx, Rx) ;                                             % Flux rate to growth (YX)
J_epg = Hill(a_epg,K_epg,Rg) ;                                             % Flux rate to product (YG)
J_epx = Hill(a_epx,K_epx,Rx) ;                                             % Flux rate to product (YX)
J_mg  = Hill(a_mg, K_mg, Rg);                                               % Flux rate to maintenance (YX)
J_mx  = Hill(a_mx, K_mx, Rx);                                               % Flux rate to maintenance (YX)
J_ecg = Hill(a_ecg,K_ecg,E) /(1+k_ecg*G);                                  % Ethanol reassimilation rate (YG)
J_ecx = Hill(a_ecx,K_ecx,E) /(1+k_ecx*X);                                  % Ethanol reassimilation rate (YX)
J_og = Hill(a_o,K_o,O);                                                    % Oxygen consumption rate (YG)
J_ox = Hill(a_o,K_o,O);                                                    % Oxygen consumption rate (YX)
J_bc = Hill(a_bc,K_bc,Cx);                                                 % Basal NAD+ consumption (YX)

% Equations
dGdt = - J_g * Ng;
dXdt = - J_x * Nx;
dEdt = b_epx * J_epx * Nx + b_epg * J_epg * Ng - J_ecx * Nx - J_ecg * Ng;
dNgdt  = g_g * J_ng * Ng;
dNxdt  = g_x * J_nx * Nx;

if Ng == 0                                                                 % Keep intracellular concentrations unchanged until the
    dRgdt = 0;                                                             % strain is added.
else
    dRgdt = b_g * J_g + b_ecg * J_ecg - J_ng - J_epg - J_mg - g_g * J_ng * Rg;
end

if Nx == 0                                                                 % Keep intracellular concentrations unchanged until the
    dRxdt = 0;                                                             % strain is added.
    dCxdt = 0;
else
    dRxdt = b_x * J_x + b_ecx * J_ecx - J_nx - J_epx - J_mx - g_x * J_nx * Rx;  % note: here g_x * J_nx is the growth rate
    dCxdt = b_o * J_ox - k_xc * J_x  - J_bc - g_x * J_nx * Cx;
end

dOdt = v_b - J_ox * Nx - J_og * Ng;

dydt = [dGdt;dXdt;dEdt;dNxdt;dNgdt;dRxdt;dRgdt;dOdt;dCxdt];

end

function y = Hill(a,K,x)
y = a*x/(x+K);
end