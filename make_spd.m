function theta = make_spd(ttheta)
% original precision matrix
myeps = 0.1;
N = size(ttheta,1);
while true
    theta = myeps*diag(ones(N,1)) + ttheta;
    ev = eig(theta);
    isspd = ~sum(ev < 0.1);
    if isspd 
        break;
    end
    myeps = myeps + 0.1;
end