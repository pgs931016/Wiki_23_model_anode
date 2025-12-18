clc; clear; close all

addpath("D:\250421연구실컴퓨터\fig_pol\RPT_OCV");
folder_path = 'D:\250421연구실컴퓨터\fig_pol\RPT_OCV';
file_list = dir(fullfile(folder_path, '*Merged.mat'));
nFiles = length(file_list);

colors = [
    0.125490196078431 0.521568627450980 0.305882352941177
    0.301960784313725 0.733333333333333 0.835294117647059
    0.572549019607843 0.368627450980392 0.623529411764706
    0.803921568627451 0.325490196078431 0.298039215686275
];

I_1C = 55.6;
kB   = 8.617e-5;

%% 저장셀 지정
para_global = cell(nFiles,1);
rmse_global = zeros(nFiles,1);
file_names  = cell(nFiles,1);

figure(1); hold on;

for f = 1:nFiles


    color_now = colors(mod(f-1,4)+1 , :);
    filePath = fullfile(folder_path, file_list(f).name);
    disp(['Loading: ', file_list(f).name]);
    data = load(filePath);
    file_names{f} = file_list(f).name;

    find_index = find([data.data_merged.OCVflag] == 1);
    data_merged = data.data_merged(find_index);

    %% 실험 파라미터 지정
    Q_all = []; t_all_raw = [];
    for i = 1:length(data_merged)
        if isfield(data_merged(i),'type') && data_merged(i).type=='C'
            t_all_raw = [t_all_raw; data_merged(i).t(end)/3600];
            Q_all     = [Q_all;     data_merged(i).Q(end)];
        end
    end

    Q0    = Q_all(1);
    Q_exp = Q_all / Q0;

    valid = (Q_all>0) & (t_all_raw>0);
    t_exp = t_all_raw(valid);
    Q_exp = Q_exp(valid);

    [t_exp, idx] = sort(t_exp);
    Q_exp = Q_exp(idx);

    %% 구간 분할
    if numel(t_exp) < 8
        num_intervals = 1;
        idx_all = [1 numel(t_exp)+1];
    else
        num_intervals = 7;
        idx_all = round(linspace(1,numel(t_exp)+1,num_intervals+1));
    end

    Q_cell = cell(num_intervals,1);
    t_cell = cell(num_intervals,1);

    for ii = 1:num_intervals
        Q_cell{ii} = Q_exp(idx_all(ii):idx_all(ii+1)-1);
        t_cell{ii} = t_exp(idx_all(ii):idx_all(ii+1)-1);
    end

    %% 시간 누적
    t_offset = 0;
    for ii = 1:num_intervals
        tt = t_cell{ii};
        tt = tt - tt(1) + t_offset;
        t_cell{ii} = tt;
        t_offset = tt(end);
    end

    %% OCV 레퍼런스
    load("H:\공유 드라이브\BSL_Data4 (2)\HNE_agedcell_2025_processed\RPT_DCIR_processed\HNE_RPT_fresh_4_1_postprocessing_HPPC.mat");
    SOC_OCV_ref = NE_OCV_linear.SOC(:);
    V_OCV_ref   = NE_OCV_linear.V(:);

    %% 전역 피팅 수행
    I_all = []; Vp_all = []; eta_all = [];
    t_all = []; dQdt_all = [];

    for ii = 1:num_intervals
        if numel(Q_cell{ii}) < 2, continue; end

        t_loc = t_cell{ii};
        dQdt  = gradient(Q_cell{ii}) ./ gradient(t_loc);
        npt   = numel(t_loc);

        SOC_vec = linspace(min(SOC_OCV_ref),max(SOC_OCV_ref),npt)';
        OCV_vec = interp1(SOC_OCV_ref,V_OCV_ref,SOC_vec,'linear','extrap');

        Vp  = 0.1 * ones(npt,1);                  eta = 0.05 * ones(npt,1);
        I   = I_1C * ones(npt,1);

        I_all    = [I_all;    I];
        Vp_all   = [Vp_all;   Vp];
        eta_all  = [eta_all;  eta];
        t_all    = [t_all;    t_loc];
        dQdt_all = [dQdt_all; dQdt];
    end

    %% 시간 격자 통일
    t_grid = linspace(min(t_all), max(t_all), 3000)';

    I_g    = interp1(t_all,I_all,t_grid,'linear','extrap');
    Vp_g   = interp1(t_all,Vp_all,t_grid,'linear','extrap');
    eta_g  = interp1(t_all,eta_all,t_grid,'linear','extrap');
    dQdt_g = interp1(t_all,dQdt_all,t_grid,'linear','extrap');

    %% 전역 초기 파라미터
    para0 = [1e-4, 3, 10, 1e-4, 10];
    lb    = [0,    0,  0,  0,    0];
    ub    = [1,   10, 50, 10,  100];

    options = optimset('display','iter','MaxIter',2000,'TolFun',1e-10);
    ms = MultiStart('UseParallel',true);

    problem = createOptimProblem('fmincon',...
        'objective',@(p) func_cost(dQdt_g,I_g,Vp_g,eta_g,t_grid,p,kB),...
        'x0',para0,'lb',lb,'ub',ub,'options',options);

    [p_hat,fval] = run(ms,problem,100);

    para_global{f} = p_hat;
    rmse_global(f) = fval;

    %% 용량 누적
    dQdt_model = func_Q(I_g,Vp_g,eta_g,t_grid,p_hat);

    Q_model = Q_exp(1) + cumtrapz(t_grid, dQdt_model);
    Q_model = Q_model / Q_model(1);

    %% 플랏
    scatter(t_exp, Q_exp,'o','MarkerEdgeColor',color_now);
    plot(t_grid, Q_model,'-','Color',color_now,'LineWidth',1.5);

end

xlabel('Time (h)');
ylabel('Capacity retention');
grid on; box on;

% save('polarization_coin_anode_GLOBAL_withModel.mat',...
%      'para_global','rmse_global','file_names');

%% 함수지정
function cost = func_cost(dQdt,I,Vp,eta,t,p,kB)
    dQ_model = func_Q(I,Vp,eta,t,p);
    cost = sqrt(mean((dQdt - dQ_model).^2));
end

function dQ = func_Q(I,Vp,eta,t,p)
    alpha=p(1); beta1=p(2); beta2=p(3);
    gamma=p(4); delta=p(5);

    t_norm = (t - t(1)) / (t(end) - t(1) + eps);

    dQ = (-alpha.*exp(beta1.*Vp) ...
          - gamma.*(I.^2).*exp(beta2.*eta)) ...
          .* exp(-delta.*sqrt(t_norm));
end


