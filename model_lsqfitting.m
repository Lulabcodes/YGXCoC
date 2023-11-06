%%% Parameter_Optimization

Para_Name = ["\alpha_x", "K_x", "\alpha_{gc}", "K_{gc}", "\alpha_{gt}", "K_{gt}", "k_x", "\alpha_{epg} and \alpha_{epx}", "K_{epg} and K_{epx}",...
"\alpha_{ecg} and \alpha_{ecx}", "K_{ecg} and K_{ecx}",   "k_{ecx}", "k_{ecg}", "\alpha_{ng} and \alpha_{nx}", "K_{ng} and K_{nx}", "m_g and m_x", "\beta_x",...
"\beta_g", "\beta_{epg} and \beta_{epx}", "\beta_{ecg} and \beta_{ecx}", "\gamma_{g} and \gamma_{x}", "v_b", "\alpha_o", "K_o", "\beta_o", "k_{xc}",...
"K_{C}", "\alpha_{bc}", "K_{bc}"];

% Input data
Data_compo = EXP_Compo;
Data_tempo = EXP_Tempo;
data = [Data_compo;Data_tempo];

% Starting points for parameter optimization, Para_0
Para_0 = Para_initial;  

% Setting boundaries for optimization
LB_Para = 0.01 * Para_0;  % Parameters are generally estimated from 0.01-100 of Para_0
UB_Para = 100 * Para_0;
UB_Para(17) = 0.98; UB_Para(18) = 0.98; % Yielding constants have ranges from 0 to the theoretical maximum
UB_Para(19) = 0.52; 
UB_Para(20) = 0.96; 
LB_Para(17:20) = 0;

% Setting optimization option
opts = optimset('MaxFunEvals',1e4);

% Solving the problem and save the results
Para_set = fmincon(@(Params)objective_function(Params, data), Para_0,[],[],[],[],LB_Para,UB_Para,[],opts);
obj = objective_function(Para_set, data); % show the optimized value of objective function
%%
%%% Calculate confidence interval of parameters
 [ci, halfwidth]= Para_confidence_interval(Para_set,data);


%%% One-at-a-time sensitivity analysis
range = 5.^(-1:0.01:1);
obj_OAT = NaN * ones(length(range),length(Para_set));

for i = 1:length(Para_set)
    Para_test = Para_set;
    Para_varied = Para_set(i) * range;
    
    for j = 1:length(Para_varied)
        Para_test(i) = Para_varied(j);
        obj_OAT(j,i) = objective_function(Para_test, data);
    end
end

save('Para_estimate.mat');

%%% Plot Figures in the main text and supplementary information
load('Para_estimate.mat');
model_plot_Fig6(Para_set);
model_plot_SFig9_SFig10(Para_set, data);


%%% Plot results of sensitivity analysis
figure;
for i = 1:length(Para_set)
    subplot(5,6,i)
    semilogx(range, sqrt(obj_OAT(:,i))/sqrt(523),'k-','LineWidth',1.5);hold on;
    title(Para_Name(i));
    ylabel('RMSD');
    xlabel('Relative change')
    ylim([0.1 0.4]);
    xlim([0.2 5]);
    xticks([0.2 1 5]);
    xticklabels([0.2 1 5]);
    ax = gca;
    ax.LineWidth = 1.0;
    ax.TickLength = [0.040,0.050];
    box on;
end
% title('One-at-a-time sensitivity analysis');