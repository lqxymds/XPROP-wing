clc
clearvars
close all

%% Example PROWIM Pylon

alpha = 0;  % Try 0 or 4
J = 0.8;
Vinf = 40;

% Wing
wing = Wing('PROWIMPylon',50);

% Propeller
beta = 23; % sets the propeller angle at r = 0.75R (by default)
prop = Propeller('PROWIM',beta);
prop.y0 = 0.325; % Location is different than on full wing PROWIM model

% Wing operating state
wingStateName = ['V = ',num2str(Vinf),', alpha = ',num2str(alpha)];
wingState = WingState(wing,wingStateName,Vinf,alpha);

% Prop operating state
propState = PropState(prop,['PROWIM - J = ',num2str(J)],J,60*Vinf/(J*prop.D));

% Slipstream
orientation = 'inboard up'; % choose the orientation of the propeller
slipstream = Slipstream('PROWIM','J','0.8'); % Import slipstream
slipstream.y0 = prop.y0; % set the location of the slipstream
slipstream = slipstream.setRadius(prop.R); % set the radius of the slipstream
slipstream = slipstream.setOrientation(orientation); % Set the orientation of the slipstream
slipstream.xLE = -0.2018; % needed for VLM Alba
slipstream.state = propState; % needed for VLM Alba
slipstream.prop = prop; % needed for VLM Alba

% Numerical results
numLL_propon = PropWingResults('LL in slipstream',wing,wingState,slipstream);
numVLM_propon = PropWingResults('Tornado',wing,wingState,slipstream);
numVLMAlba_propon = PropWingResults('VLM Alba',wing,wingState,slipstream);

% Experimental results
exp_propon = PropWingResults('Import',wing,wingState,num2str(J),num2str(alpha));
exp_propon.type = 'Experiment';
exp_propon.prop = prop;

%% Plot

figure
hold all
exp_propon.pltcompare('wing','yn','cl',...
    numLL_propon,...
    numVLM_propon,...
    numVLMAlba_propon);
exp_propon.pltfinish
