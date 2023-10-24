% Bayesian inference through MCMC simulation was used to evaluate the
% parameter uncertainty in our study.


addpath("mcmcstat-master\"); % Require the toolbox from https://github.com/mjlaine/mcmcstat
addpath("..\");

load("..\Para_estimate.mat","Para_set"); % load the model parameters, which are estimate from the least-square method.
load("..\Para_estimate.mat","data"); % load the data.
%%

%%% Calculate the mean-square-error, which is used as the initial sigma2 for MCMC simulation.
ymodel = NaN*ones(150,4);
ydata = NaN*ones(150,4);

c = 0;
for i  = 1:length(data)
    ymodel(c+1:c+length(data{i}.ydata),:) = Dolmodel(data{i},Para_set);
    ydata(c+1:c+length(data{i}.ydata),:) = data{i}.ydata(:,2:end);
    c = c+length(data{i}.ydata);
end

mse = sum((ymodel(:,[2,1,3,4])-ydata).^2)/length(ymodel);

%%% Initial set for MCMC simulation
% Parameter name
Para_Name = ["\alpha_x", "K_x", "\alpha_{gc}", "K_{gc}", "\alpha_{gt}", "K_{gt}", "k_x", "\alpha_{epg} and \alpha_{epx}", "K_{epg} and K_{epx}",...
"\alpha_{ecg} and \alpha_{ecx}", "K_{ecg} and K_{ecx}",   "k_{ecx}", "k_{ecg}", "\alpha_{ng} and \alpha_{nx}", "K_{ng} and K_{nx}", "m_g and m_x", "\beta_x",...
"\beta_g", "\beta_{epg} and \beta_{epx}", "\beta_{ecg} and \beta_{ecx}", "\gamma_{g} and \gamma_{x}", "v_b", "\alpha_o", "K_o", "\beta_o", "k_{xc}",...
"K_{C}", "\alpha_{bc}", "K_{bc}"];

% Parameter boundaries
Para_LB = zeros(1,29);
Para_UB = ones(1,29)*Inf;
Para_LB(7) = 0.01; Para_LB(12:13) = 5;
Para_UB(17:18) = 0.98; Para_UB(19) = 0.52; Para_UB(20) = 0.96;

% MCMC simulation setting
num = 2e5;
burnin = 4e4;
options.method      = 'dram';        
options.nsimu       =  num; 
options.burnintime  =  burnin;
options.adaptint    =  500; 
options.updatesigma = 1;
    
model.ssfun   = @ssDolmodel;
model.sigma2     = mse; % noninformative prior

seed = 100:100:600; 
for k = 1:length(seed)
    rng(seed(k))  % set seed for random number generator
    Para_test = Para_LB + (Para_set - Para_LB).*gamrnd(100,1/100,1,29);
    if any(Para_test>Para_UB)
        Para_test(Para_test>Para_UB) = Para_UB(Para_test>Para_UB);
    end
    
    params = cell(29,1);
    for i = 1:29
        params(i) = {{Para_Name(i), Para_test(i), Para_LB(i), Para_UB(i), NaN, inf, 1, 0}};
    end
      
    [results,chain,s2chain]=mcmcrun(model,data,params,options);

    % for each chain, sample 2000 parameter sets for evalulating model
    % prediction
    nsample = 2000;
    isample = ceil(rand(nsample,1)*(num-burnin)+burnin);
    Para_sample = chain(isample,:);
    
    pred_sample = NaN * zeros(nsample,202);

    for i = 1:nsample
        pred_sample(i,:) = model_prediction(Para_sample(i,:));
    end
    save(['mcmc_seed',num2str(seed(k))]);
end

%%% Combine multiple chains
seed = 100:100:600;
chain_total = zeros((num-burnin)*length(seed),29);
s2chain_total = zeros((num-burnin)*length(seed),4);
palchain = cell(length(seed),1);
pred_allsample = zeros(nsample*length(seed),202);

for k = 1:length(seed)
    load(['mcmc_seed',num2str(seed(k))],'chain','s2chain','pred_sample')
    chain_total((k-1)*(num-burnin)+1:k*(num-burnin),:) = chain(burnin+1:end,:);
    s2chain_total((k-1)*(num-burnin)+1:k*(num-burnin),:) = s2chain(burnin+1:end,:);
    pred_allsample((k-1)*nsample+1:k*nsample,:) = pred_sample;
    palchain{k} = chain;
end
out = mcmcpred(results,chain_total,s2chain_total,data,@Dolmodel,10000); 
save('mcmc_chain.mat', 'out', 'chain_total', 'pred_allsample','palchain', 'burnin', 'num');

%%% Display posterior distribution
LB_CI = NaN * zeros(1,29);  % The lower bound and upper bound of credible interval
UB_CI = NaN * zeros(1,29);
figure;
for i = 1: 29
    %Determine credible interval
    data_array = sort(chain_total(:,i));
    UB_CI(i) = data_array(floor(0.975*length(data_array)));   
    LB_CI(i) = data_array(floor(0.025*length(data_array)+1));
    
    % Plot figures
    subplot(6,5,i)    
    h = histogram(log10(data_array),100,'Normalization','pdf','FaceColor',"#EDB120",'FaceAlpha',0.3,'EdgeColor','none'); hold on;
    height = max(h.Values);
    plot(log10(Para_set(i)),0,'ro','MarkerFaceColor','r');
    plot(log10([LB_CI(i),LB_CI(i)]),[0 height],'r');
    plot(log10([UB_CI(i),UB_CI(i)]),[0 height],'r');
    title(strjoin([Para_Name(i),' [', num2str(log10(UB_CI(i)/LB_CI(i)),'%.2f'),']']));
    ylim([0 height])
    ax = gca;
    ax.XLim = [floor(min(ax.XLim)) floor(max(ax.XLim)+1)];
    xticks(min(ax.XLim):1+floor((max(ax.XLim)-min(ax.XLim))/6):max(ax.XLim));
    xtickformat('1e%+2.0f' );    
end

% Output a table including parameter values and confidence interval bounds.
Para_list = array2table([Para_set; LB_CI; UB_CI], 'VariableNames',Para_Name, 'RowNames',{'Value', 'Lower Bound', 'Upper Bound'});


%%% Show the correlation matrix between parameters
figure;
corr(log(chain_total));
imagesc(corr(log(chain_total(:,[1:26,28,27,29]))),[-1 1]);
colorbar;
xticks(1:29);
yticks(1:29);
xticklabels(Para_Name([1:26,28,27,29]));
yticklabels(Para_Name([1:26,28,27,29]));
%%% Calculate Gelman-Rubin-Brooks potential scale reduction factor
PSRF = psrf_plot(palchain, Para_Name, 0.05, burnin+1, num);

%%% Credible envelope for time trajectory
time_trajectory_plot(Para_set, out, data);
%%% Credible envelope for model prediction of maximal ethanol
optimal_condition_plot(Para_set, pred_allsample, 0.05)
