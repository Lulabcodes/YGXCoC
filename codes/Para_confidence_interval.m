function [ci, halfwidth]= Para_confidence_interval(Para_set,data)
% Calculate confidence interval of parameters

p = length(Para_set); % number of parameters

% number of model output = 552, instead of 600. The t=0 points are excluded as they are fixed by initial conditions.

[~, obs_residual] = modeloutput(Para_set,data);

Jm = NaN * ones(length(obs_residual), p);

step = 1e-6;
for i = 1:p
    Para_test1 = Para_set;
    Para_test1(i) = (1-step)*Para_set(i);
    Para_test2 = Para_set;
    Para_test2(i) = (1+step)*Para_set(i);    
    [out1, ~] = modeloutput(Para_test1,data);
    [out2, ~] = modeloutput(Para_test2,data);
    Jm(:,i) = (out2 - out1)/(2*step*Para_set(i));
end
ci = nlparci(Para_set,obs_residual,"Jacobian",Jm);
halfwidth = (ci(:,2)-ci(:,1))/2;
end

%%
function [modelop, obs_residual] = modeloutput(Para_set,data)

ymodel = NaN*ones(138,4);
ydata = NaN*ones(138,4);

c = 0;
for i = 1: length(data)
    Y = Dolmodel(data{i},Para_set);
    ymodel(c+1:c+length(data{i}.ydata)-1,:) = Y(2:end,:);
    ydata(c+1:c+length(data{i}.ydata)-1,:) = data{i}.ydata(2:end,[3,2,4,5]);
    c = c+length(data{i}.ydata)-1;
end
modelop = ymodel(:);
obs_residual = abs(ymodel(:) - ydata(:));
end
