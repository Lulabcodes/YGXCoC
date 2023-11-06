function Para_set = Para_initial(~)
% Initialize the parameter for least-square fitting

% assuming the only difference between specialists is substrate uptake
% pathway.

Para_set = NaN * ones(1,29);

%%% parameters for compositional DOL
a_x = 0.3;       K_x = 5;          % Xylose uptake of YX strain
Para_set(1) = a_x; Para_set(2) = K_x; 

a_g_c = 0.3;        K_g_c = 20; % Glucose uptake of YG strain
Para_set(3) = a_g_c; Para_set(4) = K_g_c;

a_g_t = 1.5;        K_g_t = 10; % Glucose uptake of YG strain
Para_set(5) = a_g_t; Para_set(6) = K_g_t;

%%% Parameter specific for xylose specialist
k_x = 0.04;         
Para_set(7) = k_x;

%%% Shared parameters
a_ep = 1;        K_ep = 0.1;  % Ethanol production
Para_set(8) = a_ep; Para_set(9) = K_ep;

a_ec = 0.1;      K_ec = 10;           k_ecx = 40;  k_ecg = 40; % Ethanol uptake
Para_set(10) = a_ec;Para_set(11) = K_ec;  Para_set(12) = k_ecx; Para_set(13) = k_ecg;

a_n = 0.01;          K_n = 0.01; % Cell growth
Para_set(14) = a_n; Para_set(15) = K_n; 

m = 0.01;     % Maintenance cost
Para_set(16) = m;       

% yielding constant
b_x = 0.60;           b_g = 0.97; 
Para_set(17) = b_x;   Para_set(18) = b_g;

b_ep = 0.51;         
Para_set(19) = b_ep; 

b_ec = 0.95;        
Para_set(20) = b_ec; 

gamma = 1;           
Para_set(21) = gamma; 

v_b = 0.01;         a_o = 0.1;        K_o = 0.1; % Oxygen-related
Para_set(22) = v_b;  Para_set(23) = a_o; Para_set(24) = K_o;

b_o = 100;          k_xc = 1;
Para_set(25) = b_o; Para_set(26) = k_xc;

K_C = 1;
Para_set(27) = K_C;

a_bc = 0.1;
Para_set(28) = a_bc;

K_bc = 0.01;
Para_set(29) = K_bc;
