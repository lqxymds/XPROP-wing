clc
clearvars
close all

%% Example PROWIM

alpha = 0; % Try 0, 4 or 8
Vinf = 50;
J = 0.85;
Tc = 0.168;

% Wing
wing = Wing('PROWIM',50);

% Propeller
beta = 25; % sets the propeller angle at r = 0.75R (by default)
prop = Propeller('PROWIM',beta);

% Prop operating state
propState = PropState(prop,['PROWIM - J = ',num2str(J)],J,60*Vinf/(J*prop.D));

wingStateName = ['V = ',num2str(Vinf),', alpha = ',num2str(alpha)];
wingState = WingState(wing,wingStateName,Vinf,alpha);

% Slipstream
orientation = 'inboard up';     % choose the orientation of the propeller
slipstream = Slipstream('PROWIM','J','0.85 Alba'); % Import slipstream
slipstream.y0 = prop.y0; % set the location of the slipstream
slipstream = slipstream.setRadius(prop.R); % set the radius of the slipstream
slipstream = slipstream.setOrientation(orientation); % Set the orientation of the slipstream


slipstream.xLE = -0.2018; % needed for VLM Alba
slipstream.state = propState; % needed for VLM Alba
slipstream.prop = prop; % needed for VLM Alba


% Numerical results
LLresults = PropWingResults('LL in slipstream',wing,wingState,slipstream);
VLMresults = PropWingResults('Tornado',wing,wingState,slipstream);
VLMAlbaresults = PropWingResults('VLM Alba',wing,wingState,slipstream);

% Experimental results
exp = PropWingResults('Import',wing,wingState,num2str(Tc),num2str(alpha));
exp.prop = prop;    % used in plotting the edges of the slipstream
exp.type = 'Experiment';

% Without propeller
woProp = PropWingResults('LL in slipstream',wing,wingState,Slipstream('empty'));
woProp.type = 'Prop off (Lifting Line)';

%% Plot
figure
hold all
exp.pltcompare('wing','yn','cl',...
    LLresults,...
    VLMresults,...
    VLMAlbaresults,...
    woProp)
exp.pltfinish()
%% Plot slipstream
% figure
% hold all
% vary = 'cdind';
% woProp.plt('wing','yn',vary);
% LLresults.plt('wing','yn',vary);
% VLMresults.plt('wing','yn',vary);
% VLMAlbaresults.plt('wing','yn',vary)
% woProp.pltfinish





