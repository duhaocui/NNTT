clear
clc
close all

f = nefHandleFunction(@func_BTT_dyn, [4, 0, 4, 0], 'isAdditive', 1);
h = nefHandleFunction(@func_rang_bear_meas, [4, 0, 2, 0], 'isAdditive', 1);