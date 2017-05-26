classdef PropResults
%   This file is part of PROP-WING TOOLS, tools for modelling the 
%   aerodynamics of propeller wing interaction.
%   This is distributed under the GNU General Public License v3 or higher.
%   Copyright (C) 2016  Kitso Epema
%
%   PropResults
%
%   Create a PropResults with
%
%   result = PropResults(prop) for an empty PropResults
%
%   result = PropResults(prop,'Import',Tagname,Tagvalue) for an imported
%               propeller result. This method searches for the results in
%               the folder:
%               PROP-WING TOOLS/Propeller/Examples
%               The filename of the results should be like:
%               Get[Propeller name without spaces]Performance.m
%               In this file the Tagname and Tagvalue are used to get the
%               correct operating conditions to import the results for.
%               It should be structured as a function. See the examples for
%               the exact details.
% 
%   result = PropResults(prop,'Import',Tagname,Tagvalue,Source)
%               Same as above, but then Source is also used to specify from
%               what results should be imported.
%
%   result = PropResults(prop,'Xrotor',state,xrotor_type) runs Xrotor based
%               on the given propeller, propState and xrotor_type.
%               Xrotor options, xrotor_type = 
%               'Graded Momentum'
%               'Potential'
%               'Vortex'
%   
%   result = PropResults(prop,'BEM Veldhuis',state) runs BEM code of
%               Veldhuis on the given propeller and propState.
%
%   Additional functionality:
%       - plt           to plot
%       - pltcompare    for a comparission between different result sources

    
    properties
        propName    % Propeller name
        prop        % Propeller
        stateName   % State name
        state       % State of the propeller
        type        % Type of analysis
        rpm
        J           % Advance ratio
        eff         % Efficiency
        Ct          % Thrust coefficient
        Cp          % Power coefficient
        Cq          % Torque coefficient
        r           % 
        va          % 
        vt          %
        aa          % Axial velocity increase factor (Vinf + va = Vinf (1 + aa))
        at          % Tangential velocity increase (Omega r - vt = Omega r (1 - at))
        vtVinf
        vaVinf
        alpha       %
        source      % source (in case of Import)
        note        % note (used in case of Import)
        beta        % Twist angle at r/R = 0.75 (used in case of Import)
        
        %TODO remove variables that are not used
        %TODO make r and rR as they should intuitively be (r dimensional,
        %       rR non-dimensional)
        %TODO rename state into propState
    end
    
    methods
        function this = PropResults(prop,type,varargin)
            %TODO, change order of input to (type,prop,varargin) to make it
            %consistent with other classes in toolbox

            if nargin < 2
                type = 'empty';
            end
            if nargin < 1
                prop = Propeller('APC11x8');
            end
            
            this.propName = prop.Name;
            this.prop = prop;
            this.type = type;
            
            switch type
                case 'empty'
                    
                case 'Import'
                    % TODO: Move ImportPerformance from propeller to
                    %           PropResults
                    if nargin > 2 
                        Tagname = varargin{1};
                        Tagvalue = varargin{2};
                        if length(varargin) > 2
                            Source = varargin{3};
                            result = prop.ImportPerformance(Tagname,Tagvalue,Source);
                        else
                            result = prop.ImportPerformance(Tagname,Tagvalue);
                        end
                    else
                        result = prop.ImportPerformance();
                    end
                    
                    this = this.saveResults(result);
                    
                case 'Actuator Disk'
                    state = varargin{1};
                    Input = varargin{2};
                    Value = varargin{3};
                    this.stateName = state.Name;
                    this.state = state;
                    result = ActuatorDisk(prop,state,Input,Value);
                    this = this.saveResults(result);
                case 'BEM'
                    state = varargin{1};
                    this.stateName = state.Name;
                    this.state = state;
                    result = BladeElement(prop,state);   
                    this = this.saveResults(result);
                    this.beta = prop.B;
                case 'Xrotor'
                    state = varargin{1};
                    this.stateName = state.Name;
                    this.state = state;
                    this.beta = prop.B;
                    if length(varargin) > 1
                        xrotor_type = varargin{2};
                    else
                        xrotor_type = 'Potential';
                    end
                    this.type = [this.type,' - ',xrotor_type];
                    result = Xrotor(prop,state,xrotor_type);
                    this = this.saveResults(result);
                case 'BEM Veldhuis'
                    state = varargin{1};
                    this.stateName = state.Name;
                    this.state = state;
                    this.beta = prop.B;
                    
                    xrotor_type = 'BEM';
                    result = Xrotor(prop,state,xrotor_type);
                    this = this.saveResults(result);
                otherwise
                    display(type)
                    error('This type of propeller analysis is not found.')
            end
        end
        function this = saveResults(this,result)
            resultNames = fieldnames(result);
            
            for n = 1:length(resultNames)
                try
                    field = char(resultNames(n));
                    this.(field) = result.(field);
                catch
                    display(field)
                    warning('Field in the results not used in PropResults')
                end
            end
        end
        function T = T(this)
            D = this.prop.R *2;
            n = this.state.n;
            T = this.Ct .* this.state.rho * D^4 * n^2;
        end
        function Ct_Vtang = Ct_Vtang(this)
            Vtang   = this.vtVinf .* this.state.Vinf;
            Q_Vtang = 0.5 * this.state.rho .* Vtang.^2;
            rr      = this.r.* this.prop.R;
            T_Vtang = trapz(rr,2 * pi .* rr .* Q_Vtang');
            D       = this.prop.R *2;
            n       = this.state.n;
            Ct_Vtang = T_Vtang ./ (this.state.rho * D^4 * n^2);
        end
        function Ct_Vaxial = Ct_Vaxial(this)
            Vaxial = (this.vaVinf + 1) .* this.state.Vinf;
            Q_Vaxial = 0.5 * this.state.rho .* Vaxial.^2;
            Q_Vaxial_upstream = 0.5 * this.state.rho .* this.state.Vinf.^2;
            rr = this.r.* this.prop.R;
            T_Vaxial = trapz(rr,2 * pi .* rr .* (Q_Vaxial' - Q_Vaxial_upstream));
            D = this.prop.R *2;
            n = this.state.n;
            Ct_Vaxial = T_Vaxial ./ (this.state.rho * D^4 * n^2);
        end
        function Ct_V = Ct_V(this)
            Va = (this.vaVinf + 1) .* this.state.Vinf;
            Vt = this.vtVinf .* this.state.Vinf;
            V2 = Va.^2 + Vt.^2;
            Q_V = 0.5 * this.state.rho .* V2;
            Q_V_upstream = 0.5 * this.state.rho .* this.state.Vinf.^2;
            rr = this.r.* this.prop.R;
            T = trapz(rr,2 * pi .* rr .* (Q_V' - Q_V_upstream));
            D = this.prop.R *2;
            n = this.state.n;
            Ct_V = T ./ (this.state.rho * D^4 * n^2);
        end
        function varargout = plt(this,Xname,Yname,LineSpec,MKsize,color)
            if nargin < 5
                MKsize = 4;
            end
            if nargin < 4
                LineSpec = 'o-';
            end 
            try
                X = this.(Xname);
            catch
                error('Could not find this.(Xname)')
            end
            try 
                Y = this.(Yname);                        
            catch
                error('Could not find this.(Yname)')
            end
%             displayName = [this.getLabelName(Yname),', ','Prop: ',this.propName,...
%                 ', (',this.type,')'];
            [a,b] = size(X);
            if a >1 && b>1
                for n = 1:length(X(1,:))
                   displayName = ['J = ',num2str(this.J(n))];
                   h = plot(X(:,n),Y(:,n),'DisplayName',displayName,...
                       'Linewidth',2,...
                        'Markersize',MKsize);
                end
            else

%                 displayName = ['\beta_{',num2str(this.prop.rB,'%1.2f'),'} = ',num2str(this.beta),' - ',this.type];
                displayName = [this.type];
                h = plot(X,Y,LineSpec,'DisplayName',displayName,...
                    'Linewidth',2,...
                    'Markersize',MKsize);
            end
            if exist('color','var')
                h.Color = color;
            end
            xlabel(this.getLabelName(Xname))
            grid on
            
            varargout{1} = h;
            
        end
        function[] = pltfinish(~)
           legend(gca,'show')
            set(gca,'fontsize',22,...
               'Linewidth',1.5)
            grid on   
        end
        function labelName = getLabelName(~,Xname)
            
            switch Xname
                case 'r'
                    labelName = 'Normalized radius r/R [-]';
                case 'vaVinf'
                    labelName = 'Normalized axial velocity V_a / V_{\infty} [-]';
                case 'vtVinf'
                    labelName = 'Normalized tangential velocity V_t / V_{\infty} [-]';
                case 'eff'
                    labelName = 'Efficiency [-]';
                case 'J'
                    labelName = 'Advance ratio J [-]';
                case 'Ct'
                    labelName = 'Thrust coefficient C_T [-]';
                case 'Cp'
                    labelName = 'Power coefficient C_P [-]';
                case 'Cq'
                    labelName = 'Torque coefficient C_Q [-]';
                otherwise
                    labelName = Xname;
            end
                    
        end
        function [] = pltcompare(this,Xname,Yname,varargin)
            Linespec = {'--o',...
                '-.^',...
                ':s',...
                '--d',...
                '-.*',...
                ':x'};
            if nargin > 3
                if length(varargin) > length(Linespec)
                    warning('pltcompare cannot handle more than 3 additional datasets')
                else
                    for n = 1:length(varargin)
                        result(n) = varargin{n};
                    end
                end
            end
            h = this.plt(Xname,Yname);
            h.LineWidth = 3;
            ylabel(this.getLabelName(Yname))
            color = get(h,'Color');

            for n = 1:length(varargin)
                if n < length(Linespec)
                    h = result(n).plt(Xname,Yname,char(Linespec(n)),10,color);
                    h.LineWidth = 1.5;
                else
                    display(['Dataset number ',n,' is not plotted.'])
                end
            end
        end
    end
end

