classdef RAC333 < NVProblem
% Defining the RAC problem for n=d=D=3
    properties
       forceReal = true;
    end
    methods
        function X = sampleOperators(self, rank)
        % the first 27 operators are the states
        % then 3 sets of 3 operators for the measurements
            dim = 3;
            X = cell(1, 36);
            % initialise random states
            for i = 1:27
                X{i} = ...
  qdimsum.Random.pureNormalizedDensityMatrix(dim);
            end
            % initialise random (rank-1 projective) measurements
            U = qdimsum.Random.unitary(dim);
            X{28} = [1 0 0; 0 0 0; 0 0 0];
            X{29} = [0 0 0; 0 1 0; 0 0 0];
            X{30} = [0 0 0; 0 0 0; 0 0 1];
            U = qdimsum.Random.unitary(dim);
            X{31} = U*[1 0 0; 0 0 0; 0 0 0]*U';
            X{32} = U*[0 0 0; 0 1 0; 0 0 0]*U';
            X{33} = U*[0 0 0; 0 0 0; 0 0 1]*U';       
            U = qdimsum.Random.unitary(dim);
            X{34} = U*[1 0 0; 0 0 0; 0 0 0]*U';
            X{35} = U*[0 0 0; 0 1 0; 0 0 0]*U';
            X{36} = U*[0 0 0; 0 0 0; 0 0 1]*U';
        end
        function K = sampleStateKraus(self)
            K = eye(3); % dimension is 3
        end
        function types = operatorTypes(self)
            types = {1:27 28:36}; % type 1: states, type 2: measurements
        end
        % action of the generators of the symmetry group on the list
        % of operators
        function generators = symmetryGroupGenerators(self)
            swapX1X2 = [1 2 3 10 11 12 19 20 21 ...
                        4 5 6 13 14 15 22 23 24 ...
                        7 8 9 16 17 18 25 26 27 ...
                        31 32 33 28 29 30 34 35 36];
            swapX1X3 = [1 10 19 4 13 22 7 16 25 ...
                        2 11 20 5 14 23 8 17 26 ...
                        3 12 21 6 15 24 9 18 27 ...
                        34 35 36 31 32 33 28 29 30];
            swapX2X3 = [1 4 7 2 5 8 3 6 9 ...
                        10 13 16 11 14 17 12 15 18 ...
                        19 22 25 20 23 26 21 24 27 ...
                        28 29 30 34 35 36 31 32 33];
            X1swap01 = [10 11 12 13 14 15 16 17 18 ...
                        1 2 3 4 5 6 7 8 9 ...
                        19 20 21 22 23 24 25 26 27 ...
                        29 28 30 31 32 33 34 35 36];
            X1swap02 = [19 20 21 22 23 24 25 26 27 ...
                        10 11 12 13 14 15 16 17 18 ...
                        1 2 3 4 5 6 7 8 9 ...
                        30 29 28 31 32 33 34 35 36];
            X1swap12 = [1 2 3 4 5 6 7 8 9 ...
                        19 20 21 22 23 24 25 26 27 ...
                        10 11 12 13 14 15 16 17 18 ...
                        28 30 29 31 32 33 34 35 36];
             generators = [swapX1X2
                           swapX1X3
                           swapX2X3
                           X1swap01
                           X1swap02
                           X1swap12];
        end
        % objective function (average success probability)
        function obj = computeObjective(self, X, K)
            obj = 0;
            for x1 = 1:3
                for x2 = 1:3
                    for x3 = 1:3
                        for y = 1:3
                          if y == 1
                             b = x1;
                         elseif y == 2
                             b = x2;
                          else
                             b = x3;
                          end
                          rho = X{x3+(x2-1)*3+(x1-1)*9};
                          M = X{27+b+(y-1)*3};
                          obj = obj + real(trace(M * rho)/81);
                        end
                    end
                end
            end
        end
    end
end
