function [functions, Ms, Ns, timings2] = speed_makegraph()

load speed1 
timings = timings1;
Ms = sizes(:, 1); 
Ns = sizes(:, 2); 
make_graph(Ms, Ns, functions, timings)

load speed2 
timings = timings2;
make_graph(Ms, Ns, functions, timings)

function make_graph(Ms, Ns, functions, timings)
figure('Position', [100, 100, 1500, 300])
hold on
fN = length(functions); 
colors = {'r', 'g', 'b'}; 

% extract min and max of mean timings: 
mn_t = Inf; 
mx_t = 0; 
for k = 1:length(Ms) 
    for l=1:fN
        t = timings{k, l}; 
        mn_t = min( mean(t), mn_t);
        mx_t = max( mean(t), mx_t);     
    end
end

for k = 1:length(Ms) 
    M = Ms(k); 
    N = Ns(k); 
    h = zeros(size(functions)); 
    for l=1:fN
        t = timings{k, l}; 
        mn = min(t) ;
        mx = max(t) ;
        mean_ = mean(t); 
        
        x = k+(l-1)*.1; 
        %h(l) = errorbar(x, log10(mean_), log10(mean_)-log10(mn), log10(mx)-log10(mean_), 'color', colors{l})
        % display just mean to avoid clutter:
        h(l)=plot(x, log10(mean_), 'x', 'color', colors{l}, ...
                  'markersize', 14, 'linew', 2); 
    end

    % plot delimiters 
    plot([k+.5, k+.5], [log10(mn_t), log10(mx_t)], 'k--'); 
end
% plot delimiter before 1st item: 
plot([0+.5, 0+.5], [log10(mn_t), log10(mx_t)], 'k--'); 

legend(h, 'LSA', 'munkres', 'lapjv', 'location', 'northwest')
ax = gca; 
N = length(Ms); 
set(ax, 'xtick', 1:N); 
% generate labels: 
labels = cell(1,N); 
for k=1:N
    labels{k} = sprintf('%ix%i', Ms(k), Ns(k)); 
end
set(ax, 'xticklab', labels); 

ylabel('log10(time in seconds)')
xlim([-1, N+1])
