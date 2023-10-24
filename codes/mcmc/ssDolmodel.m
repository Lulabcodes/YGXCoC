function obj = ssDolmodel(Para_set, data)
% sum-of-square error

ymodel = NaN*ones(150,4);
ydata = NaN*ones(150,4);

c = 0;
for i  = 1:length(data)
    ymodel(c+1:c+length(data{i}.ydata),:) = Dolmodel(data{i},Para_set);
    ydata(c+1:c+length(data{i}.ydata),:) = data{i}.ydata(:,2:end);
    c = c+length(data{i}.ydata);
end

obj = sum((ymodel(:,[2,1,3,4])-ydata).^2); 

end
