PWD = pwd;

semicolon = ';';

subdirs = {'.','misc','rv','functions','systems','estimators','identification',['help' filesep 'demo_examples']};
dirs_to_add = strcat([PWD,filesep],subdirs,semicolon);
addpath(strcat(dirs_to_add{:}));

clear PWD filesep semicolon subdirs dirs_to_add
