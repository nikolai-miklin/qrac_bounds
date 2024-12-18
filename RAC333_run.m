yalmip('clear') % clear memory
settings = NVSettings;
problem = RAC333;
% defining monomials as full levels of the NPA hierarchy
monomials = {'npa' 1};
% alternatively, defining monomials as product of types
% type 1: states, type 2: measurements
%monomials = {'families' [] [1] [2] [1 2] [1 1] [1 2 2]};
format long
% run the SDP with symmetry settings 'none', 'isotypic', 'irreps' or
% 'blocks'
upperBoundSDP = nvOptimize(problem, monomials, 'blocks', settings);
fprintf('Upper bound on the ASR: %.8f \n', upperBoundSDP)