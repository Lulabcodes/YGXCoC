function optimal_condition_plot(Para_set, pred_sample, a)

pred_ref = model_prediction(Para_set); % prediction from least-square method

nsample = length(pred_sample);
llim = ceil(nsample * a/2);
ulim = ceil(nsample * (1-a/2));

ratio_range = 0:0.01:1;
ratio_delay = -50:50;

op_r = NaN * ones(nsample,1);
op_t = NaN * ones(nsample,1);
for i = 1:nsample
    pred_com = pred_sample(i,1:length(ratio_range));
    pred_tem = pred_sample(i,length(ratio_range)+1:end);
    op_r(i) = ratio_range(pred_com==max(pred_com));
    op_t(i) = ratio_delay(pred_tem==max(pred_tem));
end
pred_sort = sort(pred_sample);


%%% Joint distribution of optimal YG fraction and delay time

figure;hold on;
hist3([op_r,op_t],'Ctrs',{0:0.01:1 -50:50},'CdataMode','auto');
axis([0.3 0.6 10 25])
% figure;hold on;
h3 = hist3([op_r,op_t],'Ctrs',{0:0.01:1 -50:50},'CdataMode','auto');

% pcolor(0:0.01:1, -50:50, h3'); %shading interp; colorbar;
view(2);
axis([0.3 0.6 10 26])
plot3(0.42, 19,1000, 'o', 'MarkerFaceColor','r', 'MarkerEdgeColor','none','MarkerSize',10);
ylabel('Optimal Delay Time (h)')
xlabel('Optimal Y_G Fraction')

for bd_th = 0:max(h3,[],'all')
    if sum(h3(h3<bd_th),'all')>0.05*length(op_r)
        break
    end
end
M = contour(0:0.01:1, -50:50, h3',[bd_th bd_th]);
plot3(M(1,2:1+M(2,1)),M(2,2:1+M(2,1)),100*ones(M(2,1)),'y-','LineWidth',2);

figure;
%%% The prediction for maximal ethanol production with varied initial YG fraction
Pexpc = [
0.1   36.8	0.296810658
0.3   45.6	0.366599598
0.4   46.9	0.368131868
0.5   47.0	0.369956744
0.7   43.3	0.356232003
0.9   41.1	0.333739837
]; % The data of validation experiments. each row shows 'YG fraction' 'Maximal ethanol production' and 'Maximal ethanol yield'
RATIO_G = (0:0.01:1)';
ini_Xyl_c = 47;
ini_Glu_c = 78;

Em_com_ref = pred_ref(1:length(RATIO_G));
subplot(1,2,1)
plot(RATIO_G, Em_com_ref/(ini_Glu_c+ini_Xyl_c),'k','LineWidth',2);hold on;
plot(Pexpc(:,1),Pexpc(:,3),'o','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','none');
ylim([0.1 0.4])
ylabel('Max yield (g/g)')
xlabel('Initial Y_G Fraction');
ax = gca;
ax.LineWidth = 1.0;
ax.TickLength = [0.020,0.025];
box on;
% Plot credible envolope

Em_com = [pred_sort(llim,1:length(RATIO_G)),pred_sort(ulim,length(RATIO_G):-1:1)];
fill([RATIO_G;RATIO_G(end:-1:1)],Em_com/(ini_Glu_c+ini_Xyl_c),...
'k','EdgeColor','None','FaceAlpha',0.3');

%%% The prediction for maximal ethanol production with varied delay time
ini_Glu_t = 71;
ini_Xyl_t = 41;
Em_tem_ref = pred_ref(length(RATIO_G)+1:end);

Pexpt = [
-32	31.321055	0.283100738
-12	31.17411	0.286896141
-6	30.135	0.284713255
0	35.902391	0.323040378
6	44.094175	0.382663568
12	48.4034	0.421744905
18	44.235645	0.37335369
27	47.87861	0.376423565
40	36.375	0.316466959
]; % The data of validation experiments. each row shows 'delay time' 'Maximal ethanol production' and 'Maximal ethanol yield'


subplot(1,2,2); hold on;
Range_T = -50:50;
plot(Range_T, Em_tem_ref/(ini_Glu_t+ini_Xyl_t),'k','LineWidth',2);

plot(Pexpt(1:end-2,1),Pexpt(1:end-2,3),'^','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','b','MarkerEdgeColor','none');
plot(Pexpt(end-1:end,1),Pexpt(end-1:end,3),'o','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','none');
ylim([0.25 0.45]);
xlim([-50 50])
ylabel('Max yield (g/g)')
xlabel('Delay Time of Y_G(h)');

ax = gca;
ax.LineWidth = 1.0;
ax.TickLength = [0.020,0.025];
box on;

% Plot credible envolope

Em_tem = [pred_sort(llim,length(RATIO_G)+1:end),pred_sort(ulim,end:-1:length(RATIO_G)+1)];
fill([Range_T,Range_T(end:-1:1)],Em_tem/(ini_Glu_t+ini_Xyl_t),...
    'k','EdgeColor','None','FaceAlpha',0.3');
