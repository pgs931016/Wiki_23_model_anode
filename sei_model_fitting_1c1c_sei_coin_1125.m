clc;clear;close all

addpath("D:\250421연구실컴퓨터\fig7_논문코드");
% D:\250421연구실컴퓨터\fig7_논문코드
folder_path = 'D:\250421연구실컴퓨터\fig_pol\RPT_OCV';
file_list = dir(fullfile(folder_path, '*Merged.mat'));

colors = [
    0.125490196078431 0.521568627450980 0.305882352941177;  % 초록
    0.301960784313725 0.733333333333333 0.835294117647059;  % 하늘
    0.572549019607843 0.368627450980392 0.623529411764706;  % 자주
    0.803921568627451 0.325490196078431 0.298039215686275   % 빨강
];


figure(1); hold on; grid on;
for f = 1:length(file_list)

    color_f = colors(mod(f-1,4)+1 , :);   

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

do_fit = true;

if f == 2 && i == 5
    do_fit = false;
end


    Q0 = data_merged(1).Q;
    Q_exp = Q / Q0;  

    tot_t = linspace(min(t), max(t), 1000);

    plot(t, Q_exp, 'o','Color', color_f,'DisplayName', 'Experimental Data'); hold on;

if do_fit == true
   
    para0 = [0.32, 1.48e2]; 

    options = optimset('display','iter', 'MaxIter',1000, 'MaxFunEvals',1e5, ...
                       'TolFun',1e-10, 'TolX',1e-10, 'FinDiffType','central');
    ms = MultiStart('Display', 'iter', 'UseParallel', true); 

    problem = createOptimProblem('fmincon', ...
        'objective', @(para) func_cost(para, t, Q_exp, kB), ... % 
        'x0', para0, 'lb', [0 0], 'ub', para0*3, 'options', options);

    num_start_points = 1000;
    [para_hat, fval, exitflag, output, solutions] = run(ms, problem, num_start_points);
    para_opt = para_hat;

    Ea_opt = para_opt(1);
    B_opt = para_opt(2);
    fprintf('Optimized Ea: %.5f eV\n', Ea_opt);
    fprintf('Optimized B: %.5e\n', B_opt);

    T = 298.15; 
    A_opt = B_opt * exp(-Ea_opt / (kB * T));

    Q_fit = 1 - A_opt * sqrt(tot_t);
    plot(tot_t, Q_fit,'-', 'LineWidth', 1,'Color', color_f, 'DisplayName', sprintf('Fitted Curve (A=%.5f)', A_opt));
end


    ylim([0.4, 1]); 
    yticks(0.4:0.1:1);
    % xticks([23.3486666666667 323.477444444445 1012.77952777778 1278.34161111111 1964.14677777778])          
    % xticklabels({'0','400','600','800','1000'})  
    xlabel('Time (h)')
    ylabel('Capacity Retention')
    h = legend('0.5C data', '0.5C fit', '1C data', '1C fit', '2C data', '2C fit','4C data', '4C fit','Location','SouthWest');
    % h = legend('QC1C Experimental', 'QC1C Predicted','Location','SouthWest');
    set(h, ...
    'FontSize',4, ...              
    'Interpreter','none', ...      
    'ItemTokenSize',[12 5]);     

    grid on; box on;
    %title('Optimized Capacity Retention Curve_1c');

end


% save('coinsei.mat','para_opt','tot_t');
save('coinsei_para_all_251128.mat');
filename = 'coinsei_edit';
figuresettings6('coinsei_edit',1200);

%%---------



function cost = func_cost(para, t, Q_exp, kB)
Ea = para(1);
B = para(2);

T = 298.15;
A = B * exp(-Ea / (kB * T));
Q_model = 1 - A * sqrt(t); 

cost = sqrt(mean((Q_exp - Q_model).^2));
end

