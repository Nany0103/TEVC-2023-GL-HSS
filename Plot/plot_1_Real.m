clc; clear;
obj = [9,9,10];
pro = {'RE91','DDMOP1','DDMOP4'};
DPlus = [7,11,13];
alpha = 1000;
selNum = 100;
numSol2 = 5000;
runNum = 21;
ylim = {[25.608,25.620],[25.34,25.42],[25.6280,25.6290],[0.244,0.252],[0.0976,0.0988],[1.38,1.50]; ...
    [57.64,57.66],[57.34,57.44],[57.6645,57.6651],[0.077,0.081],[0.028,0.0285],[0.75,0.85]};
x = {[10,100,1000,10000,100000],[10,100,1000,10000],[10,100,1000,10000,100000],[10,100,1000,10000],[10,100,1000,10000],[10,100,1000]; ...
    [10,100,1000,10000,100000,1000000],[10,100,1000,10000,100000],[10,100,1000,10000,100000,1000000],[10,100,1000,10000],[10,100,1000,10000,100000],[10,100,1000]};
HV = zeros(runNum,3,length(pro)*2);
RT = zeros(runNum,3,length(pro)*2);
for proInd = 1:length(pro)
    M = obj(proInd);
    H = getH(selNum,M);
    r = 1+1/H;
    ref = (1+1/H)*ones(1,M);
    fileName = sprintf('./Data/Approximated PF/ApproximatedPF_%s_M%d.mat',pro{proInd},M);
    PF = load(fileName).data;
    fmax = max(PF);
    fmin = min(PF);
    rowInd = proInd;
    hv_GA100 = zeros(1,runNum); rt_GA100 = zeros(1,runNum);
    hv_TGA = zeros(1,runNum); rt_TGA = zeros(1,runNum);
    hvDL4_TGA = cell(1,runNum); rtDL4_TGA = cell(1,runNum);
    hvGLNAGP_TGA = cell(1,runNum); rtGLNAGP_TGA = cell(1,runNum);
    hvGLNAGO_TGA = cell(1,runNum); rtGLNAGO_TGA = cell(1,runNum);
    hvDL4_TGA_end = zeros(1,runNum); rtDL4_TGA_end = zeros(1,runNum);
    hvLGI = zeros(1,runNum); rtLGI = zeros(1,runNum);
    for runInd = 1:runNum
        j = 100;
        % load GA100
        [hv,rt] = load_GA(ref,M,pro{proInd},selNum,runInd,fmin,fmax);
        hv_GA100(runInd) = hv(1); rt_GA100(runInd) = rt(1);
        % load TGA+DL4: Proposed one
        [hv,rt] = load_GL_TGA(ref,M,pro{proInd},selNum,alpha,numSol2,runInd,fmin,fmax);
        hv_TGA(runInd) = hv(1); rt_TGA(runInd) = rt(1);
        hvDL4_TGA{runInd} = hv; rtDL4_TGA{runInd} = rt;
        hvDL4_TGA_end(runInd) = hv(end); rtDL4_TGA_end(runInd) = rt(end);
        % load LGI-HSS
        [hvLGI(runInd),rtLGI(runInd)] = load_LGI(ref,M,pro{proInd},selNum,runInd,fmin,fmax);
        % load TGA+GLNAGO
        [hv,rt] = load_GLNAGO_TGA(ref,M,pro{proInd},selNum,numSol2,runInd,fmin,fmax);
        hvGLNAGO_TGA{runInd} = hv; rtGLNAGO_TGA{runInd} = rt;
        % save data in table
        HV(runInd,1,rowInd) = hvGLNAGO_TGA{runInd}(end); RT(runInd,1,rowInd) = rtGLNAGO_TGA{runInd}(end);
        HV(runInd,2,rowInd) = hv_GA100(runInd); RT(runInd,2,rowInd) = rt_GA100(runInd);
        HV(runInd,3,rowInd) = hvLGI(runInd); RT(runInd,3,rowInd) = rtLGI(runInd);
        %% Plot Figure
        f = figure('Position',[200,200,650,450]);
        hold on;
        fs = 21;
        lw = 1.5; ms = 20;
        plot(rtLGI(runInd),hvLGI(runInd),'s','MarkerFaceColor',0.5*ones(1,3),'MarkerSize',ms,'MarkerEdgeColor','k','LineWidth',lw);
        plot(rt_GA100(runInd),hv_GA100(runInd),'s','MarkerFaceColor','r','MarkerSize',ms,'MarkerEdgeColor','k','LineWidth',lw);
        plot(rt_TGA(runInd),hv_TGA(runInd),'s','MarkerFaceColor','y','MarkerSize',ms,'MarkerEdgeColor','k','LineWidth',lw);
        plot(rtDL4_TGA{runInd},hvDL4_TGA{runInd},'-o','MarkerFaceColor','g','MarkerSize',ms/2,'MarkerEdgeColor','k','Color','k','LineWidth',lw);
        plot(rtGLNAGO_TGA{runInd},hvGLNAGO_TGA{runInd},'-o','MarkerFaceColor','b','MarkerSize',ms/2,'MarkerEdgeColor','k','Color','k','LineWidth',lw);
        ax = gca;
        ax.XScale = 'linear';
        ax.FontName = 'Times New Roman';
        ax.FontSize = fs;
%             ax.YLim = [hv_GA100(runInd),hvLGI(runInd)];
        xlabel('Computation time (s)');
        ylabel('Hypervolume');
        ax.XGrid = 'on';
        legend({'LGI-HSS','GAHSS','TGAHSS','GL-HSS-Old','GL-HSS'}, ...
            'FontName','Times New Roman','Location','eastoutside');
        yTick = get(gca,'yTick');
        % yTickLabel = arrayfun(@(x) sprintf(['%.',num2str(p(objInd,proInd)),'f'],x), yTick,'uniformoutput',false);
%             set(gca, 'YTickLabel', yTickLabel);
        set(gca,'LineWidth', 2);
        box off;
        ax2 = axes('Position',get(gca,'Position'),...
               'XAxisLocation','top',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','k','YColor','k');
        set(ax2,'YTick', []);
        set(ax2,'XTick', []);
        set(ax2,'LineWidth', 2);
        % output high-quality vector diagram
        set(gcf, 'renderer', 'painters');
        title([pro{proInd},' ({\it m} = ',num2str(M),') run = ',num2str(runInd)], ...
            'FontWeight','normal','Interpreter','latex','FontSize',fs, ...
            'FontName','Times New Roman');
        saveas(f,sprintf('./Figure/Overall_EMO/M%d_%s_%d.png',M,pro{proInd},runInd));
        close all;
    end
end
write_table(HV,RT,pro,obj);

function write_table(HV,RT,pro,obj)
    HVAve = zeros(size(HV,3),size(HV,2));
    RTAve = zeros(size(RT,3),size(RT,2));
    for rowInd = 1:size(HV,3)
        HVAve(rowInd,:) = mean(HV(:,:,rowInd));
        RTAve(rowInd,:) = mean(RT(:,:,rowInd));
    end
    HVBet = zeros(size(HV,3),size(HV,2));
    HVWor = zeros(size(HV,3),size(HV,2));
    RTBet = zeros(size(RT,3),size(RT,2));
    RTWor = zeros(size(RT,3),size(RT,2));
    for rowInd = 1:size(HV,3)
        for runInd = 1:size(HV,1)
            for colInd = [2,3]
                if HV(runInd,colInd,rowInd)>=HV(runInd,1,rowInd)
                    HVBet(rowInd,colInd) = HVBet(rowInd,colInd)+1;
                else
                    HVWor(rowInd,colInd) = HVWor(rowInd,colInd)+1;
                end
                if RT(runInd,colInd,rowInd)<=RT(runInd,1,rowInd)
                    RTBet(rowInd,colInd) = RTBet(rowInd,colInd)+1;
                else
                    RTWor(rowInd,colInd) = RTWor(rowInd,colInd)+1;
                end
            end
        end
    end
    FID = fopen('./HV_Overall_Comparison_EMO.tex','w');
    for rowInd = 1:size(HV,3)
        if mod(rowInd-1,2)==0
            fprintf(FID,'\\multirow{2}{*}{%s}',pro{floor((rowInd-1)/2)+1});
        else
            fprintf(FID,'{}');
        end
        fprintf(FID,'&%d',obj(mod((rowInd-1),2)+1));
        for colInd = 1:size(HV,2)
            if colInd~=1
                fprintf(FID,'&%.2e (%d/%d)',HVAve(rowInd,colInd),HVBet(rowInd,colInd),HVWor(rowInd,colInd));
            else
                fprintf(FID,'&%.2e',HVAve(rowInd,colInd));
            end
        end
        fprintf(FID,'\\\\ \n');
        if mod(rowInd-1,2)==0
            fprintf(FID,'\\cline{2-5} \n');
        else
            fprintf(FID,'\\hline \n');
        end
    end
    fclose(FID);
    FID = fopen('./RT_Overall_Comparison_EMO.tex','w');
    for rowInd = 1:size(HV,3)
        if mod(rowInd-1,2)==0
            fprintf(FID,'\\multirow{2}{*}{%s}',pro{floor((rowInd-1)/2)+1});
        else
            fprintf(FID,'{}');
        end
        fprintf(FID,'&%d',obj(mod((rowInd-1),2)+1));
        for colInd = 1:size(HV,2)
            if colInd~=1
                fprintf(FID,'&%.2e (%d/%d)',RTAve(rowInd,colInd),RTBet(rowInd,colInd),RTWor(rowInd,colInd));
            else
                fprintf(FID,'&%.2e',RTAve(rowInd,colInd));
            end
        end
        fprintf(FID,'\\\\ \n');
        if mod(rowInd-1,2)==0
            fprintf(FID,'\\cline{2-5} \n');
        else
            fprintf(FID,'\\hline \n');
        end
    end
    fclose(FID);
end

function[hv,rt] = load_GL_TGA(ref,M,proDispName,selNum,alpha,numSol2,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/TGAHSS_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat', ...
        M,proDispName,selNum,numSol2,runInd);
    rtGA = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/GLHSS_TGA_M%d_%s_selNum=%d_alpha=%d_numSol2=%d_runInd=%d.mat', ...
        M,proDispName,selNum,alpha,numSol2,runInd);
    dataSelSet = load(fileName).dataSelSet;
    rt = load(fileName).runTime;
    rt = rtGA+rt;
    fileName = sprintf('./Result/Result New/HV_GLHSS_TGA_M%d_%s_selNum=%d_alpha=%d_numSol2=%d_runInd=%d.mat', ...
        M,proDispName,selNum,alpha,numSol2,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        hv = zeros(1,length(rt));
        for i=1:length(rt)
            dataSelNorm = (dataSelSet(:,:,i)-fmin)./(fmax-fmin);
            hv(i) = stk_dominatedhv(dataSelNorm,ref);
        end
        save(fileName,'hv');
    end
end

function[hv,rt] = load_GLNAGP_TGA(ref,M,proDispName,selNum,numSol2,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/TGAHSS_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat', ...
        M,proDispName,selNum,numSol2,runInd);
    rtGA = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/GLHSSNAGP_TGA_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat', ...
        M,proDispName,selNum,numSol2,runInd);
    dataSelSet = load(fileName).dataSelSet;
    rt = load(fileName).runTime;
    rt = rtGA+rt;
    fileName = sprintf('./Result/Result New/HV_GLHSSNAGP_TGA_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat', ...
        M,proDispName,selNum,numSol2,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        hv = zeros(1,length(rt));
        for i=1:length(rt)
            dataSelNorm = (dataSelSet(:,:,i)-fmin)./(fmax-fmin);
            hv(i) = stk_dominatedhv(dataSelNorm,ref);
        end
        save(fileName,'hv');
    end
end

function[hv,rt] = load_GLNAGO_TGA(ref,M,proDispName,selNum,numSol2,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/TGAHSS_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat', ...
        M,proDispName,selNum,numSol2,runInd);
    rtGA = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/GLHSSNAGO_TGA_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat', ...
        M,proDispName,selNum,numSol2,runInd);
    dataSelSet = load(fileName).dataSelSet;
    rt = load(fileName).runTime;
    rt = rtGA+rt;
    fileName = sprintf('./Result/Result New/HV_GLHSSNAGO_TGA_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat', ...
        M,proDispName,selNum,numSol2,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        hv = zeros(1,length(rt));
        for i=1:length(rt)
            dataSelNorm = (dataSelSet(:,:,i)-fmin)./(fmax-fmin);
            hv(i) = stk_dominatedhv(dataSelNorm,ref);
        end
        save(fileName,'hv');
    end
end

function[hv,rt] = load_GA(ref,M,proDispName,selNum,runInd,fmin,fmax)
    j = 100;
    fileName = sprintf('./Result/Result New/GAHSS%d_M%d_%s_selNum=%d_runInd=%d.mat', ...
        j,M,proDispName,selNum,runInd);
    rt = load(fileName).runTime;
    dataSel = load(fileName).dataSel;
    fileName = sprintf('./Result/Result New/HV_GAHSS%d_M%d_%s_selNum=%d_runInd=%d.mat', ...
        j,M,proDispName,selNum,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        dataSelNorm = (dataSel-fmin)./(fmax-fmin);
        hv = stk_dominatedhv(dataSelNorm,ref);
        save(fileName,'hv');
    end
end
function[hv,rt] = load_LGI(ref,M,proDispName,selNum,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/LGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSel = load(fileName).dataSel;
    rt = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/HV_LGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
%     if exist(fileName,'file')
%         hv = load(fileName).hv;
%     else
        dataSelNorm = (dataSel-fmin)./(fmax-fmin);
        hv = stk_dominatedhv(dataSelNorm,ref);
        save(fileName,'hv');
%     end
end

function v = nvector(m)
    for i=1:size(m,2)
        mt = m;
        mt(:,i) = [];
        v(i) = (-1)^(i-1)*det(mt);
    end
end
function H = getH(N,M)
    H = 1;
    while nchoosek(H+M,M-1) <= N
        H = H + 1;
    end
end
function h = hvc(data,i,ref)
    s = data(i,:);
    data(i,:) = [];
    datap = max(data,s);
    h = prod(ref-datap)-stk_dominatedhv(datap,ref);
end