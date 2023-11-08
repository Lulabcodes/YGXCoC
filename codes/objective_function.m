function obj = objective_function(Para_set, data)
% sum-of-square error
obj = 0;

for i  = 1:length(data)
    ymodel = Dolmodel(data{i},Para_set);
    ydata = data{i}.ydata(:,2:5); 
    hj = [39,28,24,22];
    obj = obj + sum(((ymodel-ydata)./hj).^2,'all');
end
%obj
end