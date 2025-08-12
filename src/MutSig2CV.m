function G = MutSig_2CV_v3_11_wrapper(maf,outdir,P)

demand_file(maf);
ede(outdir);

if exist('P','var')
  if ischar(P)
    params_file = P;
    P = [];
    P = process_params_file(P, params_file);
  elseif isstruct(P)
    params_file = [outdir '/params.txt'];
    write_params_file(P,params_file);
  else
    error('unknown format for P');
  end
else
  params_file = 'none';
end

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'mutation_blacklist_file',        'none');
P = impose_default_value(P,'mutation_type_dictionary_file',  'MutSig2CV_reference/mutation_type_dictionary.v4.txt');
P = impose_default_value(P,'coverage_models_mat_file',       'MutSig2CV_reference/coverage_models.agilent.mat');
P = impose_default_value(P,'basewise_coverage_fwb_file',     'none');
P = impose_default_value(P,'target_list_file',               'MutSig2CV_reference/target_list.gencode25.txt');
P = impose_default_value(P,'context_and_effect_fwb_file',    'MutSig2CV_reference/context_and_effect_track.gencodeV39.fwb');
P = impose_default_value(P,'context_and_effect_categs_file', 'MutSig2CV_reference/context_and_effect.c65e29b.txt');
P = impose_default_value(P,'covariates_file',                'MutSig2CV_reference/covariates.agilent.txt');
P = impose_default_value(P,'conservation_fwb_file',          'MutSig2CV_reference/conservation_phyloP100way_hg38.fwb');
P = impose_default_value(P,'FixedWidthBinary_jar_file',      'MutSig2CV_reference/FixedWidthBinary.jar');

P = impose_default_value(P,'build','hg38');
P = impose_default_value(P,'impute_full_coverage', true);

args = {...
    maf ...
    outdir ...
    P.mutation_blacklist_file ...
    P.mutation_type_dictionary_file ...
    P.coverage_models_mat_file ...
    P.basewise_coverage_fwb_file ...
    P.target_list_file ...
    P.context_and_effect_fwb_file ...
    P.context_and_effect_categs_file ...
    P.covariates_file ...
    P.conservation_fwb_file ...
    P.FixedWidthBinary_jar_file ...
    P.build ...
    params_file ...
    P.impute_full_coverage ...
};

G = MutSig_2CV_v3_11_core(args{:});

end
