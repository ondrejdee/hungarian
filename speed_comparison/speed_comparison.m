function [functions, Ms, Ns, timings2] = speed_comparison()

check_existence_lapjv_munkres(); 

rand('seed', 1) % set seed for repeatable results
[functions, sizes, timings1] = speed_comparison_(30000, @random1);
save speed1 functions sizes timings1

rand('seed', 1)
[functions, sizes, timings2] = speed_comparison_(30000, @random2);
save speed2 functions sizes timings2

function cost_matrix = random1_(dim, M, N, sc, shift)
% distance between two sets of random points 
% dim = dimensionality 
% M, N: cardinality of the two sets
% sc, shift: see code
u = rand(dim, M); 
v = sc*rand(dim, N) + shift; 
D = sqrt(dist2(u, v)); 
cost_matrix = D; 

function cost_matrix = random1(M, N)
cost_matrix = random1_(2, M, N, .1, 0.5); 

function cost_matrix = random2(M, N)
cost_matrix = rand(M, N); 

function [functions, sizes, timings] = speed_comparison_(repC, data_fn)

functions = {@linear_sum_assignment, @munkres, @lapjv}

% size of problems: 
sizes = [[10, 10], 
         [20, 20], 
         [50, 50], 
         [100, 100], 
         [10, 30], 
         [20, 60], 
         [50, 150], 
         [100, 300], 
         [10, 50], 
         [20, 100], 
         [50, 250], 
         [100, 500], 
         [10, 70], 
         [20, 140], 
         [50, 350], 
         [100, 700], 
         [10, 100], 
         [20, 200], 
         [50, 500], 
         [100, 1000]]; 
         

fnN = length(functions); 
sizesN = size(sizes, 1); 

timings = cell(sizesN, fnN); 

   
for i_MN = 1:sizesN 
    M = sizes(i_MN, 1);
    N = sizes(i_MN, 2);
    fprintf('%i, %ix%i\n', i_MN, M, N ); 

    % how many repetitions? 
    % we want at least 3
    repetitionN = max(3, ceil(repC/(N*M))); 
    
    t = zeros(repetitionN, fnN); 
    
    for i_rep = 1:repetitionN 
        cost_matrix = data_fn(M, N); 
        % if the matrix is rectangular, we want it to have more 
        % cols than rows: 
        assert(size(cost_matrix, 1) <= size(cost_matrix, 2));
        
        % for all functions: 
        for i_fn = 1:fnN
            fn = functions{i_fn}; 
            
            tic
            fn(cost_matrix);
            t(i_rep, i_fn) = toc; 
        end
        
        
    end
    
    % keep the timings for this setting: 
    for i_fn = 1:fnN
        timings{i_MN, i_fn} = t(:, i_fn); 
    end
end







function d2=dist2(x,y)
% function d2=dist2(x,y)
%
% INPUTS: 
% x: D-by-M matrix which stacks M vectors of dimension D
% y: D-by-N matrix which stacks N vectors of dimension D
%
% OUTPUT: 
% d2: M-by-N matrix
%     element (k,l) is a squared distance between 
%     k-th vector in x and l-th vector in y

xN=size(x,2); 
yN=size(y,2); 

d2=repmat(sum(x.^2),[yN 1])' + repmat(sum(y.^2),[xN 1]) - 2*x'*y;


