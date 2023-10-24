function PSRF = psrf_plot(palchain, Para_Name, alpha, s, e)

[PSRF,~,~]=psrf(palchain, alpha, s, e);
[PSRF, index] = sort(PSRF);
figure; bar(PSRF-1);
xticks(1:29);
xticklabels(Para_Name(index));
ylim([-1 1.0]);
yticks(-1:0.5:1.0);
yticklabels([0:0.5:2.0]);

end

