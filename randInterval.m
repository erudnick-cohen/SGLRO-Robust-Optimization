function out = randInterval(lb, ub)
    out = rand(size(lb)).*(ub-lb)+lb;
end