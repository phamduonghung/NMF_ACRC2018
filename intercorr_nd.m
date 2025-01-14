%% Intercorrelation between Hmod1 and i*Hmod2 possible after nmf

function [intercorr, timeS] = intercorr_nd(Hmod1, Hmod2, Fe)

% Intercorr normalized
% Time in seconds

[intercorr, time] = xcorr(Hmod1, Hmod2, 'unbiased');
intercorr=intercorr/(std(Hmod1)*std(Hmod2));
timeS=time/Fe;