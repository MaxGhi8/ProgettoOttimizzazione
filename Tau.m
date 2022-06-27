function [tau] = Tau(m,M)
tau = M*ones(ceil(m/M), 1);
resto = sum(tau) - m;
for i = 1:resto
    tau(i) = tau(i) - 1;
end
end