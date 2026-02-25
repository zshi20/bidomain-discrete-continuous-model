function S = sample_from_probability(P)
S = rand(size(P))<P;
end