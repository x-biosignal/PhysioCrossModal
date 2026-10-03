args = argv();
if numel(args) ~= 3
  error(["usage: generate_nonparametric_granger_fieldtrip_reference.m " ...
         "INPUT.mat OUTPUT.mat FIELDTRIP_DIR"]);
end

input_path = args{1};
output_path = args{2};
fieldtrip_dir = args{3};
addpath(fieldtrip_dir);

load(input_path);
[H, Z, Sf] = sfactorization_wilson(
  S, freq, 1000, 1e-8, "none", "chol", true, false
);

n_frequency = numel(freq);
H_for_gc = reshape(H, [1, 2, 2, n_frequency]);
Z_for_gc = reshape(Z, [1, 2, 2]);
S_for_gc = reshape(Sf, [1, 2, 2, n_frequency]);
granger = ft_connectivity_granger(
  H_for_gc, Z_for_gc, S_for_gc,
  "dimord", "rpt_chan_chan_freq", "method", "granger"
);

gc_xy = reshape(granger(1, 2, :), [1, n_frequency]);
gc_yx = reshape(granger(2, 1, :), [1, n_frequency]);
octave_version = version();
save("-mat7-binary", output_path, "H", "Z", "Sf", "gc_xy", "gc_yx",
     "octave_version");
