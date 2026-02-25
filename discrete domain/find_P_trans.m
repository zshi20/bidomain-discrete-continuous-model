function P_trans = find_P_trans(k_cat, cx, delta_t, c_50)

P_trans = 1 - exp(-k_cat .* cx .* delta_t./(c_50 + cx));

end