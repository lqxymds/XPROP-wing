clc
clearvars
close all

%% Test PropResults class

beta = 27.5;
prop = Propeller('N250',beta);

import = PropResults(prop,'Import','beta',num2str(beta));

J = import.J;
state = PropState(prop,'RPM 10000',J,import.state.rpm);

BEM = PropResults(prop,'BEM Veldhuis',state);

%% Plot
figure
hold all
BEM.pltcompare('J','eff',...
    import)
BEM.pltfinish