function [yhat] = func_sigmoid(beta,X)
yhat = 100./(1+10.^((beta(1)-X)));
end

