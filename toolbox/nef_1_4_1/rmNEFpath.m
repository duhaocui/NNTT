PWD = pwd;

semicolon = ';';

subdirs = {'.','misc','rv','functions','systems','estimators',['help' filesep 'demo_examples']};
dirs_to_remove = strcat([PWD,filesep],subdirs,semicolon);
rmpath(strcat(dirs_to_remove{:}));

clear PWD filesep semicolon subdirs dirs_to_remove