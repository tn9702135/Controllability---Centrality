%======================================================================
% MATLAB Script Information
%======================================================================
% Paper Title: "Exploring Relationships Between Controllability Criteria
% and Centrality Metrics in Brain Structural Network"
%
% Description:
%   Computes Average Controllability (AC) and Modal Controllability (MC)
%   of brain structural networks, performs correlations with various network
%   centrality metrics, calculates ICC, and optionally plots figures.
%
% MATLAB Version: R2019b
% Author: Tahereh Niyazmand
% Date: 2026-02-01
%======================================================================

%% ======================== Initialization ===========================
clear; clc; close all;

% Load dataset
data = load('NCTfMRI30SubScale60_ROI_volcorrected.mat'); % Dataset: Matrix_129, Lausanne atlas
Plot_AC_MC_figures = 1;
Epsilon = 0.025; % Coupling scaling
alpha = 0.1;     % Katz centrality attenuation factor
beta  = 1;       % Katz centrality weight for immediate neighbors
Threshold = 0.05; % Significance threshold for plots

% Scale adjacency matrices
A_all = Epsilon * data.X_ROI_volscaled;
[numSubjects,numNodes, numRaw] = size(A_all);

% Convert to cell array of square adjacency matrices
for i= 1: numSubjects
    for j= 1:numRaw
        for k =1: numNodes
            SC{i}(j,k)= A_all(i,j,k);
        end
    end
end

%% ==================== Controllability Computation =================
AC = zeros(numNodes, numSubjects);
MC = zeros(numNodes, numSubjects);

for k = 1:numSubjects
    MC(:,k) = modal_control(SC{k});
    AC(:,k) = ave_control(SC{k});
end

%% ======================= Centrality Metrics =======================
Degree      = zeros(numNodes, numSubjects);
Closeness   = zeros(numNodes, numSubjects);
Betweenness = zeros(numNodes, numSubjects);
Pagerank    = zeros(numNodes, numSubjects);
Eigenvector = zeros(numNodes, numSubjects);
Harmonic    = zeros(numNodes, numSubjects);
Katz        = zeros(numNodes, numSubjects);

for k = 1:numSubjects
    G = graph(SC{k});
    Degree(:,k)      = centrality(G,'degree','Importance',G.Edges.Weight);
    Closeness(:,k)   = centrality(G,'closeness');
    Betweenness(:,k) = centrality(G,'betweenness');
    Pagerank(:,k)    = centrality(G,'pagerank','Importance',G.Edges.Weight);
    Eigenvector(:,k) = centrality(G,'eigenvector','Importance',G.Edges.Weight);
    Harmonic(:,k)    = harmonic_centrality(SC{k});
    Katz(:,k)        = katz_centrality(SC{k}, alpha, beta);
end

%% ======================= Within-Subject Spearman ==================
Centralities = {'Degree','Closeness','Betweenness','Pagerank','Eigenvector','Harmonic','Katz'};
S = numSubjects;

median_AC_all = zeros(1,length(Centralities));
iqr_AC_all    = zeros(1,length(Centralities));
P_AC_all      = zeros(1,length(Centralities));

median_MC_all = zeros(1,length(Centralities));
iqr_MC_all    = zeros(1,length(Centralities));
P_MC_all      = zeros(1,length(Centralities));

results = struct();

for c = 1:length(Centralities)
    a = eval(Centralities{c});  % Matrix [nodes x subjects]
    rho_AC = zeros(1,S);
    rho_MC = zeros(1,S);
    
    for s = 1:S
        rho_AC(s) = corr(a(:,s), AC(:,s), 'Type','Spearman','Rows','complete');
        rho_MC(s) = corr(a(:,s), MC(:,s), 'Type','Spearman','Rows','complete');
    end
    
    median_AC_all(c) = median(rho_AC);
    iqr_AC_all(c)    = iqr(rho_AC);
    median_MC_all(c) = median(rho_MC);
    iqr_MC_all(c)    = iqr(rho_MC);
    
    [P_AC_all(c),~,~] = signrank(rho_AC);
    [P_MC_all(c),~,~] = signrank(rho_MC);
    
    results.(Centralities{c}).median_AC = median_AC_all(c);
    results.(Centralities{c}).iqr_AC    = iqr_AC_all(c);
    results.(Centralities{c}).P_AC      = P_AC_all(c);
    results.(Centralities{c}).median_MC = median_MC_all(c);
    results.(Centralities{c}).iqr_MC    = iqr_MC_all(c);
    results.(Centralities{c}).P_MC      = P_MC_all(c);
end

% Display Spearman table
T_Spearman = table(Centralities', median_AC_all', iqr_AC_all', P_AC_all', ...
                   median_MC_all', iqr_MC_all', P_MC_all', ...
                   'VariableNames', {'Centrality','Median_AC','IQR_AC','P_AC','Median_MC','IQR_MC','P_MC'});
disp('Within-subject, across-node Spearman correlation results:');
disp(T_Spearman);

%% =========================== ICC(3,1) =============================
ICC_AC = icc_3_1(AC);
ICC_MC = icc_3_1(MC);

fprintf('--------------------------------------------------------------------\n');
fprintf('ICC_AC: %.3f\n', ICC_AC);
fprintf('ICC_MC: %.3f\n\n', ICC_MC);

%% ====================== Partial Correlations =======================
partial_r_all = zeros(numSubjects,length(Centralities));
partial_p_all = zeros(numSubjects,length(Centralities));

for subj = 1:numSubjects
    MC_s = MC(:,subj);
    C = [Degree(:,subj),Closeness(:,subj),Betweenness(:,subj),Pagerank(:,subj), ...
         Eigenvector(:,subj),Harmonic(:,subj),Katz(:,subj)];
    for i = 1:length(Centralities)
        Ci = C(:,i);
        C_rest = C(:, setdiff(1:length(Centralities),i));
        [r,p] = partialcorr(MC_s, Ci, C_rest);
        partial_r_all(subj,i) = r;
        partial_p_all(subj,i) = p;
    end
end

mean_r = mean(partial_r_all);
sd_r   = std(partial_r_all);

T_partial = table(Centralities', mean_r', sd_r', 'VariableNames', {'Centrality','Mean_r','SD_r'});
disp('Partial correlations (MC vs centralities) across subjects:');
disp(T_partial);



%% ========================= Optional Plot ==========================
if Plot_AC_MC_figures
    CentralityMatrices = {Degree, Closeness, Betweenness, Pagerank, Eigenvector, Harmonic, Katz};
    CentralityNames = {'Degree','Closeness','Betweenness','Pagerank','Eigenvector','Harmonic','Katz'};
    
    for c = 1:length(CentralityMatrices)
        C_mat = CentralityMatrices{c};
        
        % Prepare figure
        figure('Name', ['AC and MC vs ' CentralityNames{c}], 'NumberTitle','off');
        
        % --- Scatter AC vs Centrality (Top-left) ---
        subplot(2,2,1); hold on;
        r_AC = zeros(1, numSubjects); p_AC = zeros(1, numSubjects);
        for i = 1:numSubjects
            scatter(C_mat(:,i), AC(:,i), 'filled');
            [r,p] = corrcoef(AC(:,i), C_mat(:,i));
            r_AC(i) = r(1,2); p_AC(i) = p(1,2);
        end
        ylabel('Average Controllability','FontWeight','bold','FontSize',12);
        xlabel(CentralityNames{c},'FontWeight','bold','FontSize',12);
        legend('the data is coloured by subject','FontWeight','bold','FontSize',10);
        box on; grid on; ax = gca; ax.FontSize = 10; ax.FontWeight = 'bold';
        hold off;
        
        % --- AC stats (Top-right, p-value log scale) ---
        subplot(2,2,2); hold on;
        analyze(p_AC, r_AC, Threshold); % This plots p-values using semilogy
        set(gca,'YScale','log');        % Explicitly set y-axis to log scale
        hold off;
        
        % --- Scatter MC vs Centrality (Bottom-left) ---
        subplot(2,2,3); hold on;
        r_MC = zeros(1, numSubjects); p_MC = zeros(1, numSubjects);
        for i = 1:numSubjects
            scatter(C_mat(:,i), MC(:,i), 'filled');
            [r,p] = corrcoef(MC(:,i), C_mat(:,i));
            r_MC(i) = r(1,2); p_MC(i) = p(1,2);
        end
        ylabel('Modal Controllability','FontWeight','bold','FontSize',12);
        xlabel(CentralityNames{c},'FontWeight','bold','FontSize',12);
        legend('the data is coloured by subject','FontWeight','bold','FontSize',10);
        box on; grid on; ax = gca; ax.FontSize = 10; ax.FontWeight = 'bold';
        hold off;
        
        % --- MC stats (Bottom-right, p-value log scale) ---
        subplot(2,2,4); hold on;
        analyze(p_MC, r_MC, Threshold); % This plots p-values using semilogy
        set(gca,'YScale','log');        % Explicitly set y-axis to log scale
        hold off;
    end
end
%% ============================ End of Script =======================

