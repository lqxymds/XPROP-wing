clc
clearvars
close all

%% Example PROWIM

alpha = 2; % Try 0, 4 or 8
Vinf = 29;
J = 0.6;
Tc = 0.3;

%% Wing
wing = Wing('Wingcontrol',50);
% wing operating state
wingStateName = ['V = ',num2str(Vinf),', alpha = ',num2str(alpha)];
wingState = WingState(wing,wingStateName,Vinf,alpha);


%% Propeller
beta = 30; % sets the propeller angle at r = 0.75R (by default)
prop = Propeller('Xprop',beta);
% Prop operating state
propStateName = ['Xprop - J = ',num2str(J)];
propState = PropState(prop,propStateName,J,60*Vinf/(J*prop.D));


%% Slipstream
orientation = 'inboard up';     % choose the orientation of the propeller
slipstream = Slipstream('Xprop','J',[num2str(J), ' Xrotor']); % Import slipstream
slipstream.y0 = prop.y0; % set span location of center of slipstream
slipstream = slipstream.setRadius(prop.R); % set the radius of the slipstream
slipstream = slipstream.setOrientation(orientation); % Set the orientation of the slipstream

slipstream.xLE = -0.195; % distance of start of slipstream (in flow direction) wrt LE
slipstream.state = propState; % needed for VLM Alba
slipstream.prop = prop; % needed for VLM Alba


% Numerical results
LLresults = PropWingResults('LL in slipstream',wing,wingState,slipstream);
VLMresults = PropWingResults('Tornado',wing,wingState,slipstream);
VLMAlbaresults = PropWingResults('VLM Alba',wing,wingState,slipstream);

% Experimental results
%exp = PropWingResults('Import',wing,wingState,num2str(Tc),num2str(alpha));
%exp.prop = prop;    % used in plotting the edges of the slipstream
%exp.type = 'Experiment';

% Without propeller
woProp = PropWingResults('LL in slipstream',wing,wingState,Slipstream('empty'));
woProp.type = 'Prop off (Lifting Line)';
woProp.prop = prop; % used in plotting the edges of the slipstream

%% Plot cl
%VLMAlbaresults
plotclflag = 0;
if plotclflag == 1
    figure
    hold all
    woProp.pltcompare('wing','yn','cl',...
        LLresults,...
        VLMresults,...
        VLMAlbaresults)
    woProp.pltfinish
end
%% Plot wing & propellel planform 
plotplanformflag = 0;
if plotplanformflag == 1
    figure
    hold on
    VLMAlbaresults.wing.plt('planform')
    prop.add2plt(VLMAlbaresults.wing)
end
%% Plot wing twist
plottwistflag = 0;
if plottwistflag == 1
    figure
    hold on
    %woProp.wing.plt('twist')
    LLresults.wing.plt('twist')
    %VLMresults.wing.plt('twist')
    VLMAlbaresults.wing.plt('twist')
end
%% Plot slipstream velocity
plotvelocityflag = 0;
if plotvelocityflag == 1
    figure
    hold on
    %woProp.slipstream.plt2()
    %LLresults.slipstream.plt2()
    %VLMresults.slipstream.plt2()
    VLMAlbaresults.slipstream.plt2()
end
%% Plot slipstream
plotslipstreamflag = 1;
if plotslipstreamflag == 1
    figure
    hold all
    vary = 'cl';
    woProp.plt('wing','yn',vary);
    LLresults.plt('wing','yn',vary);
    VLMresults.plt('wing','yn',vary);
    VLMAlbaresults.plt('wing','yn',vary)
    LLresults.pltfinish
end







