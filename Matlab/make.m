%%
src_files_base = {'cvx_matrix.c'; 'op_gradient.c'; 'op_bval.c'; 'op_beta.c'; 'op_slewrate.c'; 'op_moments.c'; 'op_eddy.c'; 'op_pns.c'; 'op_maxwell.c'; 'te_finder.c'};
%%

src_files = sprintf('../src/%s ' ,src_files_base{:});
src_files = strtrim(src_files);
src_files = ['mex_gropt_diff_fixdt.c' ' ' src_files];

inc_path = ['-I' '../src/'];
command = ['mex -v CFLAGS="$CFLAGS -std=c11" ' src_files ' ' inc_path];

eval(command);

%%
src_files = sprintf('../src/%s ' ,src_files_base{:});
src_files = strtrim(src_files);
src_files = ['mex_gropt_diff_fixN.c' ' ' src_files];

inc_path = ['-I' '../src/'];
command = ['mex -v CFLAGS="$CFLAGS -std=c11" ' src_files ' ' inc_path];

eval(command);
