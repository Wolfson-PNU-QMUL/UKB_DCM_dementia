function [p_sorted,p_corrected,H, labels, beta,SE] = UKB_DCM_dem_holmbonf(p, labels, beta,SE)

alpha = 0.05; 
m = numel(p);
H = zeros(1, m);

[p_sorted, inds] = sort(p);

p_corrected = [];
for i = 1:length(p)
    if p_sorted(i)<(alpha/(m+1-i))
        H(i) = 1;

    end
    p_corrected(i) = p_sorted(i)*(m+1-i)
end

H = logical(H);

if ~isempty(labels)
    labels = labels(inds)
else
    labels = inds;
end
beta = beta(inds);
SE = SE(inds);
end