function [e_max,R_opt,M_opt] = SeeSawQRAC(d,D,n,n_init)
% SEESAWQRAC algorithm finds in a seesaw manner a lower bound for the
% average success probability (ASP) of QRAC
%
%   Input:
%   D -- dimension of Alice's states and Bob's POVMs
%   d -- dimension of Alice's dits (normaly D=d)
%   n -- number of Alice's dits
%   n_iter -- number of iterations for each initial point
%   n_init -- number of initial points
%
%   Output:
%   e_max -- maximal found value of the ASP
%   R_opt -- states of Alice that give this maximum
%   M_opt -- measurements of Bob 

% set the parameters for the optimization
e_max = 0; % initial maximum of the ASP
e_tol = 10^-6; % target precision (used for stopping of the seesaw)
n_iter = 200; % maximum number of iterations of a single seesaw run
init = 0; % counter for the sucessful number of seesaw runs
e_all = zeros(1,n_init);
tic; % for estimating time
while init<=n_init-1 
    try % used to prevent an error from SDP solver to stop the whole program
        disp(['********* ',num2str(init+1),' : ',num2str(e_max), ' **********']);
        % first generate random states for Alice
        l_r = zeros(d^n,D^2-1);
        for i=1:d^n
            l_r(i,:) = VState(D,RandomState(D)); % vector representation of the states
        end
        % now optimize
        iter = 0;
        e_diff = 1;
        e_current = 0;
        while e_diff>e_tol && iter<=n_iter
            iter = iter+1;
            fprintf('%i: ',iter);
            [l_b,e] = opt_bob(d,D,n,l_r);
            fprintf('%f / ',e);
            [l_r,e] = opt_alice(d,D,n,l_b);
            fprintf('%f \n ',e);
            e_diff = e-e_current;
            e_current = e;
        end
        if e>=e_max
            l_r_opt = l_r;
            l_b_opt = l_b;
            e_max = e;
        end
        init=init+1;
        if init == 10
            disp([' ----------> Estimated time to finish: ',num2str(toc*n_init/10), ' sec'])
        end
        e_all(init) = e;
    catch
        fprintf('[Error] '); % to see the error message, comment out try - catch commands
    end
end
% retrieve the optimal states and measurements
for i=1:d^n 
    R_opt{i} = QState(D,[1/sqrt(D),l_r_opt(i,:)]); % turn back to matrix representation
end
for i=1:n
    M = [];
    for j=1:d
        M(:,:,j) = QState(D,l_b_opt((i-1)*d+j,:));
    end
    M_opt{i} = M;
end

% show the histogram of the outcomes
e_all = round(e_all*100)/100;
e_min = min(e_all)-0.1;
e_max_r = max(e_all);
n_int = 10; % divide all values of e on n_int intervals
N_values = zeros(1,n_int);
e_values = e_min:((e_max_r-e_min)/n_int):e_max_r; 
for k=1:n_int
    N_values(k) = numel(find(e_all>e_values(k) & e_all<=e_values(k+1)));
end

bar(fliplr(e_values(2:end)),fliplr(N_values));

end

%% SUBROUTINES

function [l_r,e] = opt_alice(d,D,n,l_b)
% OPT_ALICE is a subroutine for optimization over Alice's measurements

l_r = sdpvar(d^n,D^2-1);
R = QState(D,[1/sqrt(D),l_r(1,:)]);
Con = [R>=0];
for i=2:d^n
    R = QState(D,[1/sqrt(D),l_r(i,:)]);
    Con = [Con, R>=0];
end
% now compute the objective funciton
x = ones(1,n);
Obj = -real([1/sqrt(D),l_r(1,:)]*(sum(l_b(x+[0:d:d*(n-1)],:),1)).');
for i=2:d^n
    % for each state sum over Bob's measurements
    x = i2s(ones(1,n)*d,i);
    Obj = Obj-real([1/sqrt(D),l_r(i,:)]*(sum(l_b(x+[0:d:d*(n-1)],:),1)).');
end
sol = optimize(Con,Obj,sdpsettings('solver','mosek','verbose',0));
% If you do not use mosek, change 'mosek' to the name of the solver used
if sol.problem == 0
    l_r = value(l_r);
    e = -value(Obj)./(n*d^n);
else
    sol.info
    yalmiperror(sol.problem)
    l_r = [];
    e = [];
end
yalmip('clear')
end

function [l_b,e] = opt_bob(d,D,n,l_r)
% OPT_BOB is a subroutine for optimization over Bob's measurements

l_b = sdpvar(n*(d-1),D^2);
for y=1:n    
    l_b_e(y,:) = [sqrt(D),zeros(1,D^2-1)] - ...
        sum(l_b(((y-1)*(d-1)+1):y*(d-1),:),1);
end
M = QState(D,l_b_e(1,:));
Con = [M>=0];
for k=2:n
    M = QState(D,l_b_e(k,:));
    Con = [Con, M>=0];
end
for k=1:n*(d-1)
    M = QState(D,l_b(k,:));
    Con = [Con, M>=0];
end
% now compute the objective funciton
x = ones(1,n);
Obj = -real([1/sqrt(D),l_r(1,:)]*(sum(l_b(x+[0:(d-1):(d-1)*(n-1)],:),1)).');
for i=2:d^n
    % for each state sum over Bob's measurements
    x = i2s(ones(1,n)*d,i);
    f = find(x~=d);
    x = x+[0:(d-1):(d-1)*(n-1)];
    Obj = Obj-real([1/sqrt(D),l_r(i,:)]*(sum(l_b(x(f),:),1)).');
    Obj = Obj-real([1/sqrt(D),l_r(i,:)]*(sum(l_b_e(sdiff(1:n,f),:),1)).');
end
sol = optimize(Con,Obj,sdpsettings('solver','mosek','verbose',0));
% If you do not use mosek, change 'mosek' to the name of the solver used
if sol.problem == 0
    l_b = [value(l_b); value(l_b_e)];
    % now we need to change the order 
    ord = [reshape(1:n*(d-1),[],n); (n*(d-1)+1):n*d];
    l_b = l_b(ord(:),:);
    e = -value(Obj)./(n*d^n);
else
    sol.info
    yalmiperror(sol.problem)
    l_b = [];
    e = [];
end
yalmip('clear')
end


function R = RandomState(D)
% RANDOMSTATE generates a pure random state
psi = randn(D,1)+1i*randn(D,1);
R = psi*(psi')/(psi'*psi);
end


function r = QState(d,l)
% QSTATE is used to turn a vector representation of a density operator to
% the matrix representation
A = GMbasis(d);
r = A(:,:,1)*l(1);
for k = 2:prod(d)^2
    r = r + l(k)*A(:,:,k);
end
end

function l = VState(d,r)
% VSTATE is used to turn a matrix representation of a density operator to a
% truncated vector representation. The density operator is assumed to be
% normalized
l = zeros(1,prod(d)^2-1);
A = GMbasis(d);
for k = 2:prod(d)^2
    l(k-1) = trace(r*A(:,:,k));
end
l = real(l);
end

function S = GMbasis(d)
% GMBASIS generates Gell-Mann matrices of dimention d
S(:,:,1) = eye(d)./sqrt(d);
for l=2:d
    for k=1:l-1
        S(k,l,end+1)=1./sqrt(2);
        S(l,k,end)=1./sqrt(2);
    end
end
for l=2:d
    for k=1:l-1
        S(k,l,end+1)=-1i./sqrt(2);
        S(l,k,end)=1i./sqrt(2);
    end
end
for l=1:d-1
    S(1:l+1,1:l+1,end+1) = [eye(l), zeros(l,1); zeros(1,l), -l].*sqrt(1./(l*(l+1)));
end    
end

function s = i2s(D,in)
% I2S is essentially the same as ind2sub, but faster. 
%    D -- vector of dimensions.
%    For D = [2 2 .. 2] use I2V, which is much faster.
%
if any(in>prod(D))
    error('Index cannot exceed prod(D).');
end
if any(in<=0)
    error('Index must be a positive integer.');
end
D = fliplr(D);
s = zeros(numel(in),numel(D));
for k=1:numel(in)
    for n=numel(D)-1:-1:1
        r = rem(in(k)-1,prod(D(1:n)))+1;
        s(k,n+1) = (in(k)-r)/prod(D(1:n))+1;
        in(k) = in(k) - (s(k,n+1)-1)*prod(D(1:n));
    end
    s(k,1) = rem(in(k)-1,D(1))+1;
    s(k,:) = fliplr(s(k,:));
end
end

function c = sdiff(a,b)
% SDIFF does the same as setdiff for two lists a and b, but sdiff keeps
%    the order. 
%    Does not work for rows.
%
%    See also SETDIFF
%
[c,index_c] = setdiff(a,b);
[~,order_c] = sort(index_c);
c = c(order_c);
end

