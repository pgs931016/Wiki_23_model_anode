clc; clear; close all

addpath("D:\250421연구실컴퓨터\fig_pol\RPT_OCV");
folder_path = 'D:\250421연구실컴퓨터\fig_pol\RPT_OCV';
file_list = dir(fullfile(folder_path, '*Merged.mat'));
nFiles = length(file_list);

colors = [
    0.125490196078431 0.521568627450980 0.305882352941177;  % 초록
    0.301960784313725 0.733333333333333 0.835294117647059;  % 하늘
    0.572549019607843 0.368627450980392 0.623529411764706;  % 자주
    0.803921568627451 0.325490196078431 0.298039215686275   % 빨강
];


I_1C = 55.6;
kB   = 8.617e-5; % eV/K
T    = 298.15;   % K (고정)

para_cell = cell(nFiles, 1);   
file_names = cell(nFiles, 1); 
figure(1); hold on;


for f = 1:nFiles
        color_f = colors(mod(f-1,4)+1 , :);   % 파일마다 다른색상


        filePath = fullfile(folder_path, file_list(f).name);
        disp(['Loading: ', file_list(f).name]);
        data = load(filePath);

        find_index = find([data.data_merged.OCVflag] == 1);
        data_merged = data.data_merged(find_index);

          num_data = length(data_merged);
          t = zeros(1, num_data);
          Q = zeros(1, num_data);

          kB = 8.617e-5; % eV/K
          T = 298.15;

for i = 1:num_data
        t(i) = data_merged(i).t(end) / 3600; 
        Q(i) = data_merged(i).Q; 
end



% 대응되는 Merged7 / Merged7_1
anode_merged7_files = data_merged;

anode_merged7_1_files = "G:\공유 드라이브\BSL_Data4 (1)\HNE_agedcell_2025_processed\RPT_DCIR_processed\HNE_RPT_fresh_4_1_postprocessing_HPPC.mat";


anode_data = anode_merged7_files; 

    load(anode_merged7_1_files); 

    color_now = color_f(mod(f-1, size(color_f,1)) + 1, :);


% NE_OCV_linear는 struct로 가정: SOC, V 가 벡터
SOC_OCV_ref = NE_OCV_linear.SOC(:);
V_OCV_ref   = NE_OCV_linear.V(:);




%% data_merged에서 charge 구간('C')의 누적 t[h], Q 추출

Q_all = []; t_all_raw = []; cycle_all = [];
for i = 1:length(data_merged)
    if isfield(data_merged(i),'type') && data_merged(i).type=='C'
        t_all_raw = [t_all_raw; data_merged(i).t(end)/3600];
        Q_all     = [Q_all;     data_merged(i).Q(end)];
        cycle_all = [cycle_all; data_merged(i).cycle];
    end
end
if isempty(Q_all), error('유효한 Q 데이터가 없습니다.'); end
Q0    = Q_all(1);
Q_exp = Q_all / Q0;  


valid = (Q_all>0) & (t_all_raw>0);
Q_vals = Q_all(valid);
t_vals = t_all_raw(valid);
[t_vals, idx_sort] = sort(t_vals);
Q_vals = Q_vals(idx_sort);


if numel(Q_vals) < 8
    num_intervals = 1;
    idx_all = [1, numel(Q_vals)+1];
else
    num_intervals = 7;
    idx_all = round(linspace(1, numel(Q_vals)+1, num_intervals+1));
end

Q_cell = cell(1, num_intervals);
t_cell = cell(1, num_intervals);
for ii = 1:num_intervals
    i1 = idx_all(ii);
    i2 = idx_all(ii+1)-1;
    Q_cell{ii} = Q_vals(i1:i2);
    t_cell{ii} = t_vals(i1:i2);
end


t_all = []; t_offset = 0;
for ii = 1:num_intervals
    if isempty(t_cell{ii}), continue; end
    tt = t_cell{ii}(:);
    tt = tt - tt(1) + t_offset;
    t_cell{ii} = tt;
    t_offset   = tt(end);
    t_all      = [t_all; tt];
end


V_start = 1; V_end = length(data_merged);
pick_idx = round(linspace(V_start, V_end, num_intervals));
pick_idx = max(1, min(pick_idx, length(data_merged)));

y1 = struct([]);
OCPp_cell = cell(1,num_intervals);
OCPn_cell = cell(1,num_intervals);


for ii = 1:num_intervals
    npt = numel(Q_cell{ii});
    if npt==0, y1(ii).V=[]; y1(ii).I=[]; y1(ii).OCV=[]; continue; end

    base = pick_idx(ii);
    Vmin = min(data_merged(base).V); Vmax = max(data_merged(base).V);
    Imin = min(data_merged(base).I); Imax = max(data_merged(base).I);
    if isempty(Vmin) || isempty(Vmax) || isempty(Imin) || isempty(Imax)
        Vmin=0; Vmax=0; Imin=I_1C; Imax=I_1C;
    end

    y1(ii).V   = linspace(Vmin, Vmax, npt).';
    y1(ii).I   = linspace(Imin, Imax, npt).';


    SOC_vec = linspace(min(SOC_OCV_ref), max(SOC_OCV_ref), npt).';
    OCV_vec = interp1(SOC_OCV_ref, V_OCV_ref, SOC_vec, 'linear', 'extrap');
    y1(ii).OCV = flipud(OCV_vec);   % 원 의도 유지


    OCPp_cell{ii} = linspace(min(anode_data(ii).OCPp), max(anode_data(ii).OCPp), npt).';
    OCPn_cell{ii} = linspace(min(anode_data(ii).OCPn), max(anode_data(ii).OCPn), npt).';
end


Vp_cell  = cell(1,num_intervals);
eta_cell = cell(1,num_intervals);
I_cell   = cell(1,num_intervals);

for ii = 1:num_intervals
    npt = numel(Q_cell{ii});
    if npt==0, Vp_cell{ii}=[]; eta_cell{ii}=[]; I_cell{ii}=[]; continue; end
    Vp_cell{ii}  = y1(ii).V(:) - y1(ii).OCV(:);
    eta_cell{ii} = (OCPp_cell{ii} - (y1(ii).V(:) - OCPn_cell{ii})) - OCPn_cell{ii};
    I_cell{ii}   = y1(ii).I(:);
end


dQdt_cell = cell(1,num_intervals);
for ii = 1:num_intervals
    if numel(Q_cell{ii})<2 || numel(t_cell{ii})<2
        dQdt_cell{ii} = zeros(size(Q_cell{ii}));
    else
        dQdt_cell{ii} = gradient(Q_cell{ii}(:)) ./ gradient(t_cell{ii}(:));
    end
end


para0 = [1e-4, 3, 10, 1e-4, 10];
lb    = [0,    0,  0,  0,    0];
ub    = [1,   10, 50, 10,  100];



scatter(t_all,  Q_exp, 'o', 'MarkerEdgeColor', color_now); 


do_fit = true;
if f == 2 && i == 5
    do_fit = false;
end

if do_fit == true

for ii = 1:num_intervals
    if isempty(Q_cell{ii}) || numel(Q_cell{ii})<2
        dQdt_model_cell{ii} = []; para_cell{ii} = []; cap_model_cell{ii} = [];
        continue
    end


    t_local = t_cell{ii}(:);
    t_grid  = linspace(t_local(1), t_local(end), 1000).';

    % 모든 입력 동일 격자로 보간
    I_data    = interp1(t_local, I_cell{ii}(:),    t_grid, 'linear','extrap');
    Vp_data   = interp1(t_local, Vp_cell{ii}(:),   t_grid, 'linear','extrap');
    eta_data  = interp1(t_local, eta_cell{ii}(:),  t_grid, 'linear','extrap');
    dQdt_data = interp1(t_local, dQdt_cell{ii}(:), t_grid, 'linear','extrap');

    options = optimset('display','iter', 'MaxIter',1000, 'MaxFunEvals',1e5, ...
                       'TolFun',1e-10, 'TolX',1e-10, 'FinDiffType','central');
    ms = MultiStart('Display', 'iter', 'UseParallel', true); 

    problem = createOptimProblem('fmincon', ...
        'objective', @(p) func_cost(dQdt_data, I_data, Vp_data, eta_data, t_grid, p, kB), ...
        'x0', para0, 'lb', lb, 'ub', ub, 'options', options);

    [x_hat, fval, exitflag, output, solutions] = run(ms, problem, 100); 

    % 모델 dQ/dt
    dQdt_model = func_Q(I_data, Vp_data, eta_data, t_grid, x_hat);
    dQdt_model_cell{ii} = dQdt_model;
    para_cell{ii}       = x_hat;
    file_names{f} = file_list(f).name;
    para_cell{f} = para_cell;

    % 용량 모델 누적(구간 기준)
    if ii==1
        cap_model_cell{ii} = Q_cell{ii}(1) + cumtrapz(t_grid, dQdt_model);
    else
        cap_model_cell{ii} = cap_model_cell{ii-1}(end) + cumtrapz(t_grid, dQdt_model);
    end

    % 로컬 RMSE 
    Q_actual_i = interp1(t_local, Q_cell{ii}(:), t_grid, 'linear','extrap');
    rmse_local(ii) = sqrt(mean( (Q_actual_i - cap_model_cell{ii}(:)).^2 ));
end


Q_initial = Q_cell{1}(1);
Q_cap_all = []; t_all_concat = [];
for ii = 1:num_intervals
    if isempty(Q_cell{ii}), continue; end
    Q_cap_all    = [Q_cap_all;    Q_cell{ii}(:) ./ Q_initial];
    t_all_concat = [t_all_concat; t_cell{ii}(:)];
end

% 모델: 구간별 cap_model_cell과 해당 t_grid를 이어붙임
cap_model_concat = []; t_model_concat = [];
for ii = 1:num_intervals
    if isempty(cap_model_cell{ii}), continue; end
    ngrid  = numel(cap_model_cell{ii});
    t_loc  = t_cell{ii}(:);
    t_grid = linspace(t_loc(1), t_loc(end), ngrid).';
    cap_model_concat = [cap_model_concat; cap_model_cell{ii}(:)];
    t_model_concat   = [t_model_concat;   t_grid];
end
% 정규화(용량유지율)
cap_model_norm = cap_model_concat / cap_model_concat(1);

% 최종 x축: tot_t(1000점)
tot_t = linspace(min(t_all_concat), max(t_all_concat), 1000).';

% 실험/모델 모두 tot_t로 보간 (같은 길이!)
Q_cap_1000 = interp1(t_all_concat,  Q_cap_all,    tot_t, 'linear','extrap');
model_1000 = interp1(t_model_concat, cap_model_norm, tot_t, 'linear','extrap');

%% ===== [7] RMSE(같은 격자에서) 및 플롯 =====
Q_cap_1000 = Q_cap_1000(:);
model_1000 = model_1000(:);

rmse_abs = sqrt(mean( (Q_cap_1000 - model_1000).^2 )); % fraction
rmse_pct = rmse_abs * 100;



plot(tot_t, model_1000, '-', 'Color', color_now, 'LineWidth', 1);
end
xlabel('Time (h)');
ylabel('Capacity Retention');

 h = legend('0.5C data', '0.5C fit', '1C data', '1C fit', '2C data', '2C fit','4C data', '4C fit','Location','Southwest');

    ylim([0.4, 1]); 
    yticks(0.4:0.1:1); 

    set(h, ...
    'FontSize',4, ...              
    'Interpreter','none', ...       
    'ItemTokenSize',[12 5]);    
    grid on; box on;

% title(sprintf('RMSE = %.3g (%.2f%%)', rmse_abs, rmse_pct));
end


 figuresettings6('coin_anodepotential_251128_2',1200);
 % %figuresettings6('1C1C_anodepotential2',1200);
 save('polarization_coin_anode_251128_2.mat','para_cell',"rmse_pct",'model_1000','tot_t');
 save('anodepotential_coin_paraall_251128_2')
 %save('polarization_1c1c_anode2.mat','para_cell',"rmse_pct",'model_1000','tot_t');


function [cost] = func_cost(dQdt, I, Vp, eta_n, t, para, kB) %#ok<INUSD>
    % 입력 길이/형상 정리
    [I, Vp, eta_n, t, dQdt] = force_same_length(I, Vp, eta_n, t, dQdt);
    dQdt_model = func_Q(I, Vp, eta_n, t, para);
    w = ones(size(dQdt));
    cost = sqrt(mean( w .* (dQdt - dQdt_model).^2 ));
end

function dQdt_model = func_Q(I, Vp, eta_n, t, para)
    [I, Vp, eta_n, t] = force_same_length(I, Vp, eta_n, t);

    alpha = para(1);
    beta1 = para(2);
    beta2 = para(3);
    gamma = para(4);
    delta = para(5);

    % 시간 정규화(0~1)
    if max(t)==min(t)
        t_norm = zeros(size(t));
    else
        t_norm = (t - t(1)) / (max(t)-t(1));
    end

    dQdt_model = (-alpha .* exp(beta1 .* Vp) ...
                  - gamma .* (I.^2) .* exp(beta2 .* eta_n)) ...
                  .* exp(-delta .* sqrt(t_norm));
end

function unilen = force_same_length(varargin)
    % 모든 입력을 열벡터로 만들고 최소 길이에 맞춰 절단
    n = inf;
    tmp = cell(1,nargin);
    for k = 1:nargin
        x = varargin{k};
        x = x(:);
        tmp{k} = x;
        n = min(n, numel(x));
    end
    for k = 1:nargin
        unilen{k} = tmp{k}(1:n);
    end
end

function [SOC_OCV_ref, V_OCV_ref] = load_soc_v_from_file(matPath)
    S = load(matPath);
    fn = fieldnames(S);

    % 1) struct 중 'SOC'와 'V' 보유 탐색
    for k = 1:numel(fn)
        val = S.(fn{k});
        if isstruct(val) && all(isfield(val, {'SOC','V'}))
            SOC = val.SOC; V = val.V;
            if isnumeric(SOC) && isnumeric(V)
                SOC_OCV_ref = SOC(:); V_OCV_ref = V(:);
                return
            end
        end
    end

   
    if isfield(S,'SOC') && isfield(S,'V') && isnumeric(S.SOC) && isnumeric(S.V)
        SOC_OCV_ref = S.SOC(:); V_OCV_ref = S.V(:);
        return
    end

    patSOC = ["SOC","soc","Soc","SOC_ref","soc_ref"];
    patV   = ["V","v","OCV","ocv","Voc","Vref","OCV_ref"];
    candSOC = []; candV = [];
    for k = 1:numel(fn)
        val = S.(fn{k});
        if isnumeric(val) && isvector(val)
            name = fn{k};
            if any(strcmpi(name, patSOC)) || contains(lower(name),'soc')
                candSOC = [candSOC, k]; %#ok<AGROW>
            end
            if any(strcmpi(name, patV)) || contains(lower(name),{'ocv','voc'})
                candV = [candV, k]; %#ok<AGROW>
            end
        end
    end
    if ~isempty(candSOC) && ~isempty(candV)
        best = [];
        bestDiff = inf;
        for i = candSOC
            for j = candV
                n1 = numel(S.(fn{i})(:)); n2 = numel(S.(fn{j})(:));
                if abs(n1-n2) < bestDiff
                    bestDiff = abs(n1-n2); best = [i j];
                end
            end
        end
        if ~isempty(best)
            SOC_OCV_ref = S.(fn{best(1)})(:);
            V_OCV_ref   = S.(fn{best(2)})(:);
            return
        end
    end


end
