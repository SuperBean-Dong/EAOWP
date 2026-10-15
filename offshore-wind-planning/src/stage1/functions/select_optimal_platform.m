function [optimal_platform, min_cost, detail] = select_optimal_platform(params,life_farm,annual_interest_rate)
%SELECT_OPTIMAL_PLATFORM Site-adaptive support-platform lifecycle-cost selection.
% MATLAB implementation aligned with the final Appendix-B parameter table.
%
% Important retained modelling choices requested by the user:
%   * 30 < H <= 100 m still allows Fixed + SSP + TLP (code takes precedence).
%   * SSP mooring-line length keeps the source-code relation 450/135*H.
%
% Lifecycle present cost:
%   C_life = C_dev + C_eq + C_ins + PV(C_O&M) + PV(C_dec)
%
% Floating support-body cost uses the final MCF formulation:
%   C_mat   = m_p * c_steel
%   C_manuf = MCF_p * C_mat
%   C_sup   = n_T * (C_mat + C_manuf)
%
% O&M:
%   annual O&M = P_cap,kW * (k_p + k_d*D_port) * USD2CNY
% Decommissioning:
%   PV(C_dec) = C_ins * k_dec / (1+r)^L

validateattributes(life_farm,{'numeric'},{'scalar','positive','finite'});
validateattributes(annual_interest_rate,{'numeric'},{'scalar','nonnegative','finite'});

% Depth-dependent admissible set. Keep the current-code treatment at 30-100 m.
if params.depth <= 30
    platform_types = {'固定式'};
elseif params.depth <= 100
    platform_types = {'固定式','半潜式','张力腿'};
else
    platform_types = {'半潜式','单柱式','张力腿'};
end

costs = zeros(1,numel(platform_types));
components = cell(1,numel(platform_types));
for i = 1:numel(platform_types)
    pnow = params;
    pnow.platform_type = platform_types{i};
    [costs(i),components{i}] = platform_lifecycle_cost(pnow,life_farm,annual_interest_rate);
end

[min_cost,idx] = min(costs);
optimal_platform = platform_types{idx};

detail = struct();
detail.candidate_platforms = platform_types;
detail.candidate_cost_CNY = costs;
detail.candidate_components = components;
detail.selected_index = idx;
detail.selected_platform = optimal_platform;
detail.minimum_lifecycle_present_cost_CNY = min_cost;
detail.selected_components = components{idx};
end

function [total_cost,D] = platform_lifecycle_cost(params,L,r)
%PLATFORM_LIFECYCLE_COST Present-value lifecycle cost for one platform concept.

required = {'n_T','P_rated','distance_port','depth','turbine_unit_cost','steel_price', ...
    'GBP2CNY','USD2CNY','EUR2CNY','lb','lc','dt'};
for k = 1:numel(required)
    if ~isfield(params,required{k})
        error('select_optimal_platform:MissingParameter','Missing params.%s.',required{k});
    end
end

nT = double(params.n_T);
PrMW = double(params.P_rated);
Dport = double(params.distance_port); % km
H = double(params.depth);             % m
Pcap_kW = nT*PrMW*1000;
Pcap_MW = nT*PrMW;

%% 1. Development and licensing
if Pcap_MW < 300
    fixed_cost_GBP = 3684373;
else
    fixed_cost_GBP = 23378000;
end
survey_cost_GBP_per_MW = 97792;
C_dev = (fixed_cost_GBP + survey_cost_GBP_per_MW*Pcap_MW)*params.GBP2CNY;

%% 2. Equipment acquisition
C_turbine = Pcap_kW*params.turbine_unit_cost;

isFixed = strcmp(params.platform_type,'固定式');

% Defaults retained for reporting.
platform_mass_t = NaN;
MCF = NaN;
C_mat_unit = 0;
C_manuf_unit = 0;
C_support = 0;
C_moor = 0;
mooring_length_m = 0;
n_mooring_lines = 0;
anchor_mass_t = 0;

if ~isFixed
    switch params.platform_type
        case '半潜式' % SSP
            platform_mass_t = 696;
            MCF = 2.00;
            mooring_length_m = (450/135)*H;
            n_mooring_lines = 6;
            anchor_mass_t = 3.15;
            mooring_unit_cost_CNY_per_m = 22;
            % Final document geometry (kept explicit for traceability).
            platform_freeboard_m = 12;
            platform_draft_m = 10;
            inferior_pontoon_diameter_m = 2;
            platform_horizontal_dimension_m = 76;
        case '单柱式' % Spar
            platform_mass_t = 989;
            MCF = 1.20;
            mooring_length_m = (390/135)*H; % 2.8889H
            n_mooring_lines = 3;
            anchor_mass_t = 3.15;
            mooring_unit_cost_CNY_per_m = 1000/60;
            platform_freeboard_m = 12;
            platform_draft_m = 120;
            inferior_pontoon_diameter_m = NaN;
            platform_horizontal_dimension_m = 9.40;
        case '张力腿' % TLP
            platform_mass_t = 965;
            MCF = 1.50;
            mooring_length_m = (85/135)*H; % 0.6296H
            n_mooring_lines = 8;
            anchor_mass_t = 8.10;
            mooring_unit_cost_CNY_per_m = 1000/60;
            platform_freeboard_m = 12;
            platform_draft_m = 45;
            inferior_pontoon_diameter_m = 3;
            platform_horizontal_dimension_m = 54;
        otherwise
            error('Unknown platform type: %s',params.platform_type);
    end

    % Final Appendix-B manufacturing formulation.
    C_mat_unit = platform_mass_t*params.steel_price;
    C_manuf_unit = MCF*C_mat_unit;
    C_support = nT*(C_mat_unit + C_manuf_unit);

    % Mooring system.
    C_mooring_line = nT*n_mooring_lines*mooring_length_m*mooring_unit_cost_CNY_per_m;
    anchor_unit_cost_CNY_per_kg = 2*params.EUR2CNY; % 2 EUR/kg
    C_anchor = nT*n_mooring_lines*anchor_mass_t*1000*anchor_unit_cost_CNY_per_kg;
    C_moor = C_mooring_line + C_anchor;
else
    platform_freeboard_m = NaN;
    platform_draft_m = NaN;
    inferior_pontoon_diameter_m = NaN;
    platform_horizontal_dimension_m = NaN;
    % Depth-dependent fixed-foundation cost [EUR/kW] converted to CNY.
    C_support = (1.276*H^2 - 6.975*H + 534.2)/10000 * 472 * params.EUR2CNY * Pcap_kW;
end

C_eq = C_turbine + C_support + C_moor;

%% 3. Installation
C_port = 0;
C_transport = 0;
C_platform_install = 0;
C_moor_install = 0;

if ~isFixed
    % 3.1 Port handling.
    t_lift_turbine_h = 3;
    t_platform_port_h = 3; % final Table B4: common platform port time
    port_rent_CNY_per_m2_h = 0.02*params.EUR2CNY/24; % 0.02 EUR/m2/day
    crane_CNY_per_h = 833*params.EUR2CNY;

    % Final formula uses maximum blade chord lc rather than rotor diameter.
    A_turbine_total_m2 = 3*nT*(params.lb*params.lc + pi*(params.dt/2)^2);
    C_port_turbine = A_turbine_total_m2*t_lift_turbine_h*port_rent_CNY_per_m2_h + ...
        nT*t_lift_turbine_h*crane_CNY_per_h;

    switch params.platform_type
        case '半潜式'
            A_platform_unit_m2 = (platform_horizontal_dimension_m^2/2)*sqrt(3);
        case '单柱式'
            A_platform_unit_m2 = (platform_draft_m+platform_freeboard_m)*platform_horizontal_dimension_m;
        case '张力腿'
            A_platform_unit_m2 = (platform_draft_m+platform_freeboard_m+inferior_pontoon_diameter_m)*platform_horizontal_dimension_m;
    end
    A_platform_total_m2 = nT*A_platform_unit_m2;
    C_port_platform = A_platform_total_m2*t_platform_port_h*port_rent_CNY_per_m2_h + ...
        nT*t_platform_port_h*crane_CNY_per_h;
    C_port = C_port_turbine + C_port_platform;

    % 3.2 Marine transport: final Table B5.
    c_tug_CNY_per_h = 22500*params.EUR2CNY/24;
    c_barge_CNY_per_h = 35000*params.EUR2CNY/24;
    c_fc_without_CNY_per_h = 116000*params.EUR2CNY/24;
    c_fc_with_CNY_per_h = 811900*params.EUR2CNY/24;
    v_tug_kmh = 3.6*3.6;
    v_barge_kmh = 3.6*3.6;
    v_fc_kmh = 3.09*3.6;
    weather_factor = 0.75;
    n_turbine_per_tug = 1;
    n_turbine_per_fc = 1;

    switch params.platform_type
        case '半潜式'
            Nv = ceil(nT/n_turbine_per_tug);
            t_trans_h = 2*Dport/v_tug_kmh;
            C_transport = Nv*c_tug_CNY_per_h*t_trans_h/weather_factor;
        case '单柱式'
            Nv = ceil(nT/n_turbine_per_tug);
            t_trans_h = n_turbine_per_tug*t_lift_turbine_h + 2*Dport/v_barge_kmh;
            % Table B5 specifies tug + two barges.
            c_v_CNY_per_h = c_tug_CNY_per_h + 2*c_barge_CNY_per_h;
            C_transport = Nv*c_v_CNY_per_h*t_trans_h/weather_factor;
        case '张力腿'
            Nv = ceil(nT/n_turbine_per_fc);
            t_trans_h = n_turbine_per_fc*t_lift_turbine_h + 2*Dport/v_fc_kmh;
            C_transport = Nv*c_fc_with_CNY_per_h*t_trans_h/weather_factor;
    end

    % 3.3 Offshore platform installation (Table B6).
    % The final parameter document does not redefine both sub-time terms;
    % retain the source-code values while applying the final vessel rates.
    t_IT_h = 3;
    t_imp_h = 8;
    switch params.platform_type
        case '半潜式'
            C_platform_install = 0;
        case '单柱式'
            t_inst_h = nT*(t_IT_h+t_imp_h);
            C_platform_install = c_fc_without_CNY_per_h*t_inst_h/weather_factor;
        case '张力腿'
            t_inst_h = nT*(t_IT_h+t_imp_h);
            C_platform_install = c_fc_with_CNY_per_h*t_inst_h/weather_factor;
    end

    % 3.4 Complete mooring-and-anchor installation: 3.4 day/line.
    t_moor_per_line_h = 3.4*24;
    C_AHV_CNY_per_h = 36466*params.EUR2CNY/24;
    C_personnel_CNY_per_h = 4780;
    n_mooring_points = nT*n_mooring_lines;
    C_moor_install = n_mooring_points*(C_AHV_CNY_per_h + C_personnel_CNY_per_h)*t_moor_per_line_h;

    C_ins = C_port + C_transport + C_platform_install + C_moor_install;
else
    % Fixed-bottom transport and installation empirical expression.
    C_ins = ((2 + 2*Dport/480)*(329476/5))/10000 * params.EUR2CNY * Pcap_kW;
end

%% 4. O&M present value: final platform-specific absolute coefficients
if isFixed
    k_p = 19.7;  % USD/(kW*yr)
    k_d = 0.013; % USD/(kW*yr*km)
    k_dec = 0.70;
else
    k_p = 20.8;  % USD/(kW*yr)
    k_d = 0.022; % USD/(kW*yr*km)
    k_dec = 0.60;
end
C_om_annual = Pcap_kW*(k_p + k_d*Dport)*params.USD2CNY;
PVIFA = annuityPVFactor(L,r);
C_om_PV = C_om_annual*PVIFA;

%% 5. Decommissioning present value
C_dec_final = C_ins*k_dec;
C_dec_PV = C_dec_final/(1+r)^L;

%% Total lifecycle present cost
total_cost = C_dev + C_eq + C_ins + C_om_PV + C_dec_PV;

D = struct();
D.platform_type = params.platform_type;
D.total_capacity_MW = Pcap_MW;
D.development_licensing_CNY = C_dev;
D.turbine_procurement_CNY = C_turbine;
D.support_structure_CNY = C_support;
D.mooring_system_CNY = C_moor;
D.equipment_acquisition_CNY = C_eq;
D.port_handling_CNY = C_port;
D.marine_transport_CNY = C_transport;
D.offshore_platform_installation_CNY = C_platform_install;
D.mooring_installation_CNY = C_moor_install;
D.installation_CNY = C_ins;
D.annual_OM_CNY = C_om_annual;
D.PV_OM_CNY = C_om_PV;
D.decommissioning_final_CNY = C_dec_final;
D.PV_decommissioning_CNY = C_dec_PV;
D.lifecycle_present_cost_CNY = total_cost;
D.platform_mass_t = platform_mass_t;
D.manufacturing_complexity_factor = MCF;
D.material_cost_per_platform_CNY = C_mat_unit;
D.manufacturing_cost_per_platform_CNY = C_manuf_unit;
D.mooring_length_m = mooring_length_m;
D.n_mooring_lines = n_mooring_lines;
D.anchor_mass_t = anchor_mass_t;
D.platform_freeboard_m = platform_freeboard_m;
D.platform_draft_m = platform_draft_m;
D.inferior_pontoon_diameter_m = inferior_pontoon_diameter_m;
D.platform_horizontal_dimension_m = platform_horizontal_dimension_m;
D.OM_base_USD_per_kWyr = k_p;
D.OM_distance_USD_per_kWyrkm = k_d;
D.decommissioning_factor = k_dec;
end

function f = annuityPVFactor(L,r)
if r == 0
    f = L;
else
    f = (1-(1+r)^(-L))/r;
end
end
