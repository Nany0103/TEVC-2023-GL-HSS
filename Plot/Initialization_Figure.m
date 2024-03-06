clc; clear;
obj = [8,10];
alpha = 1000;
pro = {'DTLZ1','DTLZ2','MinusDTLZ1','MinusDTLZ2','WFG3','DTLZ7'};
proO = {@DTLZ1,@DTLZ2,@MinusDTLZ1,@MinusDTLZ2,@WFG3,@DTLZ7};
selNum = 100;
numSol2 = 5000;
p = [3,2,4,3,4,2;3,2,4,3,4,2];
runNum = 21;
ylim = {[25.608,25.620],[25.34,25.42],[25.6280,25.6290],[0.244,0.252],[0.0976,0.0988],[1.38,1.50]; ...
    [57.64,57.66],[57.34,57.44],[57.6645,57.6651],[0.077,0.081],[0.028,0.0285],[0.75,0.85]};
x = {[10,100,1000,10000,100000],[10,100,1000,10000],[10,100,1000,10000,100000],[10,100,1000,10000],[10,100,1000,10000],[10,100,1000]; ...
    [10,100,1000,10000,100000,1000000],[10,100,1000,10000,100000],[10,100,1000,10000,100000,1000000],[10,100,1000,10000],[10,100,1000,10000,100000],[10,100,1000]};
algNum = 4;
% GL-HSS, GL-HSS^IDSS, GL-HSS^DSS, GL-HSS^Rand
HVEMO = zeros(runNum,algNum,length(pro)*2);
RTEMO = zeros(runNum,algNum,length(pro)*2);
for objInd = 1:2
    M = obj(objInd);
    H = getH(selNum,M);
    r = 1+1/H;
    ref = (1+1/H)*ones(1,M);
    for proInd = 1:6
        if proInd~=5
            Pro = proO{proInd};
            Problem = Pro('M',M);
            fmax = max(Problem.optimum);
            fmin = min(Problem.optimum);
        else
            fileName = sprintf('./Data/Approximated PF/ApproximatedPF_%s_M%d.mat',pro{proInd},M);
            PF = load(fileName).data;
            fmax = max(PF);
            fmin = min(PF);
        end
        ThreeInd = (proInd-1)*length(obj)+objInd;
        for runInd = 1:runNum
            j = 100;
            % load TGA+GLNAGO
            [hv,rt] = load_GLNAGO_TGA(ref,M,pro{proInd},selNum,numSol2,runInd,fmin,fmax);
            HVEMO(runInd,1,ThreeInd) = hv(end); RTEMO(runInd,1,ThreeInd) = rt(end);
            % load Rand+GLNAGO
            [hv,rt] = load_GLNAGO_Rand(ref,M,pro{proInd},selNum,runInd,fmin,fmax);
            HVEMO(runInd,2,ThreeInd) = hv(end); RTEMO(runInd,2,ThreeInd) = rt(end);
            % load IDSS+GLNAGO
            [hv,rt] = load_GLNAGO_IDSS(ref,M,pro{proInd},selNum,runInd,fmin,fmax);
            HVEMO(runInd,3,ThreeInd) = hv(end); RTEMO(runInd,3,ThreeInd) = rt(end);
            % load DSS+GLNAGO
            [hv,rt] = load_GLNAGO_DSS(ref,M,pro{proInd},selNum,runInd,fmin,fmax);
            HVEMO(runInd,4,ThreeInd) = hv(end); RTEMO(runInd,4,ThreeInd) = rt(end);
        end
    end
end

%% load PF data
obj = [8,10];
alpha = 1000;
proEMO =  {'DTLZ1','DTLZ2','MinusDTLZ1','MinusDTLZ2','WFG3','DTLZ7'};
proEMODisp = {'DTLZ1','DTLZ2','M-DTLZ1','M-DTLZ2','WFG3','DTLZ7'};
pro = {'Linear','Concave','Convex','I-Linear','I-Concave','I-Convex'};
selNum = 100;
numSol2 = 5000;
p = [1,1,2,2,3,3,4,4,4,4,1,1,3,4,3,4,3,3,3,4,1,1,1,1,1,1,1];
runNum = 21;
ylim = {[25.608,25.620],[25.34,25.42],[25.6280,25.6290],[0.244,0.252],[0.0976,0.0988],[1.38,1.50]; ...
    [57.64,57.66],[57.34,57.44],[57.6645,57.6651],[0.077,0.081],[0.028,0.0285],[0.75,0.85]};
x = {[10,100,1000,10000,100000],[10,100,1000,10000],[10,100,1000,10000,100000],[10,100,1000,10000],[10,100,1000,10000],[10,100,1000]; ...
    [10,100,1000,10000,100000,1000000],[10,100,1000,10000,100000],[10,100,1000,10000,100000,1000000],[10,100,1000,10000],[10,100,1000,10000,100000],[10,100,1000]};
HVPF = zeros(runNum,algNum,length(pro)*2);
RTPF = zeros(runNum,algNum,length(pro)*2);
for objInd = 1:2
    M = obj(objInd);
    H = getH(selNum,M);
    r = 1+1/H;
    for proInd = 1:6
        fprintf('M%d_%s\n',M,pro{proInd});
        rowInd = (proInd-1)*length(obj)+objInd;
        ThreeInd = (proInd-1)*length(obj)+objInd;
        for runInd = 1:runNum
            ref = ones(1,M)*r;
            j = 100;
            % load TGA+GLNAGO
            [hv,rt] = load_GLNAGO_TGA_PF(ref,M,pro{proInd},selNum,numSol2,runInd);
            HVPF(runInd,1,ThreeInd) = hv(end); RTPF(runInd,1,ThreeInd) = rt(end);
            % load Rand+GLNAGO
            [hv,rt] = load_GLNAGO_Rand_PF(ref,M,pro{proInd},selNum,runInd);
            HVPF(runInd,2,ThreeInd) = hv(end); RTPF(runInd,2,ThreeInd) = rt(end);
            % load IDSS+GLNAGO
            [hv,rt] = load_GLNAGO_IDSS_PF(ref,M,pro{proInd},selNum,runInd);
            HVPF(runInd,3,ThreeInd) = hv(end); RTPF(runInd,3,ThreeInd) = rt(end);
            % load DSS+GLNAGO
            [hv,rt] = load_GLNAGO_DSS_PF(ref,M,pro{proInd},selNum,runInd);
            HVPF(runInd,4,ThreeInd) = hv(end); RTPF(runInd,4,ThreeInd) = rt(end);
        end
    end
end

obj = [9,9,10];
proReal = {'RE91','DDMOP1','DDMOP4'};
DPlus = [7,11,13];
alpha = 1000;
selNum = 100;
numSol2 = 5000;
runNum = 21;
HVReal = zeros(runNum,algNum,length(proReal));
RTReal = zeros(runNum,algNum,length(proReal));
for proInd = 1:length(proReal)
    M = obj(proInd);
    H = getH(selNum,M);
    r = 1+1/H;
    ref = (1+1/H)*ones(1,M);
    fileName = sprintf('./Data/Approximated PF/ApproximatedPF_%s_M%d.mat',proReal{proInd},M);
    PF = load(fileName).data;
    fmax = max(PF);
    fmin = min(PF);
    for runInd = 1:runNum
        j = 100;
        % load TGA+GLNAGO
        [hv,rt] = load_GLNAGO_TGA(ref,M,proReal{proInd},selNum,numSol2,runInd,fmin,fmax);
        HVReal(runInd,1,proInd) = hv(end); RTReal(runInd,1,proInd) = rt(end);
        % load Rand+GLNAGO
        [hv,rt] = load_GLNAGO_Rand(ref,M,proReal{proInd},selNum,runInd,fmin,fmax);
        HVReal(runInd,2,proInd) = hv(end); RTReal(runInd,2,proInd) = rt(end);
        % load IDSS+GLNAGO
        [hv,rt] = load_GLNAGO_IDSS(ref,M,proReal{proInd},selNum,runInd,fmin,fmax);
        HVReal(runInd,3,proInd) = hv(end); RTReal(runInd,3,proInd) = rt(end);
        % load DSS+GLNAGO
        [hv,rt] = load_GLNAGO_DSS(ref,M,proReal{proInd},selNum,runInd,fmin,fmax);
        HVReal(runInd,4,proInd) = hv(end); RTReal(runInd,4,proInd) = rt(end);
    end
end

HV = zeros(runNum,algNum,27);
HV(:,:,1:12) = HVPF; HV(:,:,13:24) = HVEMO; HV(:,:,25:27) = HVReal;
RT = zeros(runNum,algNum,27); 
RT(:,:,1:12) = RTPF; RT(:,:,13:24) = RTEMO; RT(:,:,25:27) = RTReal;
pro = {'Linear','Linear','Concave','Concave','Convex','Convex', ...
    'I-Linear','I-Linear','I-Concave','I-Concave','I-Convex','I-Convex', ...
    'DTLZ1','DTLZ1','DTLZ2','DTLZ2','MinusDTLZ1','MinusDTLZ1', ...
    'MinusDTLZ2','MinusDTLZ2','WFG3','WFG3','DTLZ7','DTLZ7','RE91','DDMOP1','DDMOP4'};
proDisp = {'Linear','Linear','Concave','Concave','Convex','Convex', ...
    'I-Linear','I-Linear','I-Concave','I-Concave','I-Convex','I-Convex', ...
    'DTLZ1','DTLZ1','DTLZ2','DTLZ2','M-DTLZ1','M-DTLZ1', ...
    'M-DTLZ2','M-DTLZ2','WFG3','WFG3','DTLZ7','DTLZ7','RE9-7-1','DDMOP1','DDMOP4'};
obj = [repmat([8,10],1,12),9,9,10];
yLimVal = {[57.33,57.43],[0.0236,0.0244],[0.067,0.073],[10.05,13]};
count = 1;
for proInd = [4,10,18,25]
    M = obj(proInd);
    hvP = HV(:,1,proInd);
    medVal = median(hvP);
    runInd = find(hvP==medVal);
    runInd = runInd(1);
    if proInd<=12
        [hv,rt] = load_GLNAGO_TGA_PF(ref,M,pro{proInd},selNum,numSol2,runInd);
        HVGL = hv; RTGL = rt;
        [hv,rt] = load_GLNAGO_Rand_PF(ref,M,pro{proInd},selNum,runInd);
        HVRand = hv; RTRand = rt;
        [hv,rt] = load_GLNAGO_IDSS_PF(ref,M,pro{proInd},selNum,runInd);
        HVIDSS = hv; RTIDSS = rt;
        [hv,rt] = load_GLNAGO_DSS_PF(ref,M,pro{proInd},selNum,runInd);
        HVDSS = hv; RTDSS = rt;
    else
        [hv,rt] = load_GLNAGO_TGA(ref,M,pro{proInd},selNum,numSol2,runInd,fmin,fmax);
        HVGL = hv; RTGL = rt;
        [hv,rt] = load_GLNAGO_Rand(ref,M,pro{proInd},selNum,runInd,fmin,fmax);
        HVRand = hv; RTRand = rt;
        [hv,rt] = load_GLNAGO_IDSS(ref,M,pro{proInd},selNum,runInd,fmin,fmax);
        HVIDSS = hv; RTIDSS = rt;
        [hv,rt] = load_GLNAGO_DSS(ref,M,pro{proInd},selNum,runInd,fmin,fmax);
        HVDSS = hv; RTDSS = rt;
    end
    f = figure('Position',[200,200,550,450]);
    hold on;
    fs = 25;
    lw = 1.5; ms = 20;
    % plot GL-HS
    plot(RTGL,HVGL,'-o','MarkerFaceColor','b','MarkerSize',ms/2,'MarkerEdgeColor','k','Color','k','LineWidth',lw);
    % plot Rand-GL-HSS
    plot(RTRand,HVRand,'-o','MarkerFaceColor','k','MarkerSize',ms/2,'MarkerEdgeColor','k','Color','k','LineWidth',lw);
    % plot IDSS-GL-HSS
    plot(RTIDSS,HVIDSS,'-o','MarkerFaceColor',[238,176,175]/255,'MarkerSize',ms/2,'MarkerEdgeColor','k','Color','k','LineWidth',lw);
    % plot DSS-GL-HSS
    plot(RTDSS,HVDSS,'-o','MarkerFaceColor',[177,196,77]/255,'MarkerSize',ms/2,'MarkerEdgeColor','k','Color','k','LineWidth',lw);
    ax = gca;
    ax.YLim = yLimVal{count};
    count = count+1;
    ax.FontName = 'Times New Roman';
    ax.FontSize = fs;
    xlabel('Computation time (s)');
    ylabel('Hypervolume (HV)');
    ax.XGrid = 'on';
    set(gca,'LineWidth', 2);
    yTick = get(gca,'yTick');
    yTickLabel = arrayfun(@(x) sprintf(['%.',num2str(p(proInd)),'f'],x), yTick,'uniformoutput',false);
    set(gca, 'YTickLabel', yTickLabel);
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
    saveas(f,sprintf('./Figure/Initialization Final/%s_M%d.emf',pro{proInd},M));
    close all;
end

write_table(HV,RT,proDisp,obj);

function write_table(HV,RT,pro,obj)
    HVAve = zeros(size(HV,3),size(HV,2));
    RTAve = zeros(size(RT,3),size(RT,2));
    HVBest = false(size(HV,3),size(HV,2));
    RTBest = false(size(HV,3),size(HV,2));
    for rowInd = 1:size(HV,3)
        HVAve(rowInd,:) = mean(HV(:,:,rowInd));
        RTAve(rowInd,:) = mean(RT(:,:,rowInd));
        ind_hv_best = find(HVAve(rowInd,:)==max(HVAve(rowInd,:)));
        ind_rt_best = find(RTAve(rowInd,:)==min(RTAve(rowInd,:)));
        HVBest(rowInd,ind_hv_best) = true;
        RTBest(rowInd,ind_rt_best) = true;
    end
    HVBet = zeros(size(HV,3),size(HV,2));
    HVWor = zeros(size(HV,3),size(HV,2));
    RTBet = zeros(size(RT,3),size(RT,2));
    RTWor = zeros(size(RT,3),size(RT,2));
    for rowInd = 1:size(HV,3)
        for runInd = 1:size(HV,1)
            for colInd = 2:size(HV,2)
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
    FID = fopen('./Table/HV_Initialization_Comparison.tex','w');
    for rowInd = 1:size(HV,3)-3
        if mod(rowInd-1,2)==0
            fprintf(FID,'\\multirow{2}{*}{%s}',pro{rowInd});
        else
            fprintf(FID,'{}');
        end
        fprintf(FID,'&%d',obj(mod((rowInd-1),2)+1));
        for colInd = 1:size(HV,2)
            if HVBet(rowInd,colInd)>HVWor(rowInd,colInd), sym = '+';
            else, sym = '-'; end
            if colInd~=1
                if HVBest(rowInd,colInd)
                    fprintf(FID,'&\\cellcolor[gray]{0.9} %.2e (%s)',HVAve(rowInd,colInd),sym);
                else
                    fprintf(FID,'&%.2e (%s)',HVAve(rowInd,colInd),sym);
                end
            else
                if HVBest(rowInd,colInd)
                    fprintf(FID,'&\\cellcolor[gray]{0.9} %.2e',HVAve(rowInd,colInd));
                else
                    fprintf(FID,'&%.2e',HVAve(rowInd,colInd));
                end
            end
        end
        fprintf(FID,'\\\\ \n');
        if mod(rowInd-1,2)==0
            fprintf(FID,'\\cline{2-%d} \n',size(HV,2)+2);
        else
            fprintf(FID,'\\hline \n');
        end
    end
    objR = [9,9,10];
    proR = {'RE9-7-1','DDMOP1','DDMOP4'};
    for rowInd = size(HV,3)-2:size(HV,3)
        fprintf(FID,'%s',proR{rowInd-24});
        fprintf(FID,'&%d',objR(rowInd-24));
        for colInd = 1:size(HV,2)
            if HVBet(rowInd,colInd)>HVWor(rowInd,colInd), sym = '+';
            else, sym = '-'; end
            if colInd~=1
                if HVBest(rowInd,colInd)
                    fprintf(FID,'&c\\cellcolor[gray]{0.9} %.2e (%s)',HVAve(rowInd,colInd),sym);
                else
                    fprintf(FID,'&%.2e (%s)',HVAve(rowInd,colInd),sym);
                end
            else
                if HVBest(rowInd,colInd)
                    fprintf(FID,'&\\cellcolor[gray]{0.9} %.2e',HVAve(rowInd,colInd));
                else
                    fprintf(FID,'&%.2e',HVAve(rowInd,colInd));
                end
            end
        end
        fprintf(FID,'\\\\ \n');
        fprintf(FID,'\\hline \n');
    end
    fprintf(FID,'\\multicolumn{3}{c|}{(+/-)}');
    for colInd = 2:size(HV,2)
        fprintf(FID,'&(%d/%d)',sum(HVBet(:,colInd)),sum(HVWor(:,colInd)));
    end
    fprintf(FID,'\\\\ \n');
    fprintf(FID,'\\hline \n');
    fclose(FID);
    FID = fopen('./Table/RT_Initialization_Comparison.tex','w');
    for rowInd = 1:size(HV,3)-3
        if mod(rowInd-1,2)==0
            fprintf(FID,'\\multirow{2}{*}{%s}',pro{rowInd});
        else
            fprintf(FID,'{}');
        end
        fprintf(FID,'&%d',obj(mod((rowInd-1),2)+1));
        for colInd = 1:size(HV,2)
            if RTBet(rowInd,colInd)>RTWor(rowInd,colInd), sym = '+';
            else, sym = '-'; end
            if colInd~=1
                if RTBest(rowInd,colInd)
                    fprintf(FID,'&\\cellcolor[gray]{0.9} %.2e (%s)',RTAve(rowInd,colInd),sym);
                else
                    fprintf(FID,'&%.2e (%s)',RTAve(rowInd,colInd),sym);
                end
            else
                if RTBest(rowInd,colInd)
                    fprintf(FID,'&\\cellcolor[gray]{0.9} %.2e',RTAve(rowInd,colInd));
                else
                    fprintf(FID,'&%.2e',RTAve(rowInd,colInd));
                end
            end
        end
        fprintf(FID,'\\\\ \n');
        if mod(rowInd-1,2)==0
            fprintf(FID,'\\cline{2-%d} \n',size(HV,2)+2);
        else
            fprintf(FID,'\\hline \n');
        end
    end
    for rowInd = size(HV,3)-2:size(HV,3)
        fprintf(FID,'%s',proR{rowInd-24});
        fprintf(FID,'&%d',objR(rowInd-24));
        for colInd = 1:size(HV,2)
            if RTBet(rowInd,colInd)>RTWor(rowInd,colInd), sym = '+';
            else, sym = '-'; end
            if colInd~=1
                if RTBest(rowInd,colInd)
                    fprintf(FID,'&\\cellcolor[gray]{0.9} %.2e (%s)',RTAve(rowInd,colInd),sym);
                else
                    fprintf(FID,'&%.2e (%s)',RTAve(rowInd,colInd),sym);
                end
            else
                if RTBest(rowInd,colInd)
                    fprintf(FID,'&\\cellcolor[gray]{0.9} %.2e',RTAve(rowInd,colInd));
                else
                    fprintf(FID,'&%.2e',RTAve(rowInd,colInd));
                end
            end
        end
        fprintf(FID,'\\\\ \n');
        fprintf(FID,'\\hline \n');
    end
    fprintf(FID,'\\multicolumn{3}{c|}{(+/-)}');
    for colInd = 2:size(HV,2)
        fprintf(FID,'&(%d/%d)',sum(RTBet(:,colInd)),sum(RTWor(:,colInd)));
    end
    fprintf(FID,'\\\\ \n');
    fprintf(FID,'\\hline \n');
    fclose(FID);
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
function[hv,rt] = load_GLNAGO_TGA_PF(ref,M,proDispName,selNum,numSol2,runInd)
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
            hv(i) = stk_dominatedhv(dataSelSet(:,:,i),ref);
        end
        save(fileName,'hv');
    end
end

function[hv,rt] = load_GLNAGO_IDSS_PF(ref,M,proDispName,selNum,runInd)
    fileName = sprintf('./Result/Result New/IDSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    rtIDSS = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/GLHSSNAGO_IDSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSelSet = load(fileName).dataSelSet;
    rt = load(fileName).runTime;
    rt = rtIDSS+rt;
    fileName = sprintf('./Result/Result New/HV_GLHSSNAGO_IDSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        hv = zeros(1,length(rt));
        for i=1:length(rt)
            hv(i) = stk_dominatedhv(dataSelSet(:,:,i),ref);
        end
        save(fileName,'hv');
    end
end

function[hv,rt] = load_GLNAGO_DSS_PF(ref,M,proDispName,selNum,runInd)
    fileName = sprintf('./Result/Result New/DSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    rtDSS = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/GLHSSNAGO_DSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSelSet = load(fileName).dataSelSet;
    rt = load(fileName).runTime;
    rt = rtDSS+rt;
    fileName = sprintf('./Result/Result New/HV_GLHSSNAGO_DSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        hv = zeros(1,length(rt));
        for i=1:length(rt)
            hv(i) = stk_dominatedhv(dataSelSet(:,:,i),ref);
        end
        save(fileName,'hv');
    end
end

function[hv,rt] = load_GLNAGO_Rand_PF(ref,M,proDispName,selNum,runInd)
    fileName = sprintf('./Result/Result New/Rand_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    rtRand = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/GLHSSNAGO_Rand_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSelSet = load(fileName).dataSelSet;
    rt = load(fileName).runTime;
    rt = rtRand+rt;
    fileName = sprintf('./Result/Result New/HV_GLHSSNAGO_RandS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        hv = zeros(1,length(rt));
        for i=1:length(rt)
            hv(i) = stk_dominatedhv(dataSelSet(:,:,i),ref);
        end
        save(fileName,'hv');
    end
end

function[hv,rt] = load_KME_PF(ref,M,proDispName,selNum,runInd)
    fileName = sprintf('./Result/Result New/KME_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSel = load(fileName).dataSel;
    rt = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/HV_KME_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        hv = stk_dominatedhv(dataSel,ref);
        save(fileName,'hv');
    end
end

function[hv,rt] = load_KME(ref,M,proDispName,selNum,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/KME_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSel = load(fileName).dataSel;
    rt = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/HV_KME_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        dataSelNorm = (dataSel-fmin)./(fmax-fmin);
        hv = stk_dominatedhv(dataSelNorm,ref);
        save(fileName,'hv');
    end
end

function[hv,rt] = load_TGI_PF(ref,M,proDispName,selNum,runInd)
    fileName = sprintf('./Result/Result New/TGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSel = load(fileName).dataSel;
    rt = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/HV_TGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        hv = stk_dominatedhv(dataSel,ref);
        save(fileName,'hv');
    end
end

function[hv,rt] = load_TGI(ref,M,proDispName,selNum,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/TGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSel = load(fileName).dataSel;
    rt = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/HV_TGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        dataSelNorm = (dataSel-fmin)./(fmax-fmin);
        hv = stk_dominatedhv(dataSelNorm,ref);
        save(fileName,'hv');
    end
end

function[hv,rt] = load_LGI_PF(ref,M,proDispName,selNum,runInd)
    fileName = sprintf('./Result/Result New/LGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSel = load(fileName).dataSel;
    rt = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/HV_LGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        hv = stk_dominatedhv(dataSel,ref);
        save(fileName,'hv');
    end
end

function[hv,rt] = load_GA_PF(ref,M,proDispName,selNum,runInd)
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
        hv = stk_dominatedhv(dataSel,ref);
        save(fileName,'hv');
    end
end

function[hv,rt] = load_APL_PF(ref,M,proDispName,selNum,runInd)
    fileName = sprintf('./Result/Result New/APLHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    rt = load(fileName).runTime;
    dataSelSet = load(fileName).dataSelSet;
    fileName = sprintf('./Result/Result New/HV_APLHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        hv = zeros(1,length(rt));
        for i=1:length(rt)
            hv(i) = stk_dominatedhv(dataSelSet(:,:,i),ref);
        end
        save(fileName,'hv');
    end
end

function[hv,rt] = load_APL(ref,M,proDispName,selNum,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/APLHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    rt = load(fileName).runTime;
    dataSelSet = load(fileName).dataSelSet;
    fileName = sprintf('./Result/Result New/HV_APLHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
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

function[hv,rt] = load_GLNAGO_Rand(ref,M,proDispName,selNum,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/Rand_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    rtRand = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/GLHSSNAGO_Rand_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSelSet = load(fileName).dataSelSet;
    rt = load(fileName).runTime;
    rt = rtRand+rt;
    fileName = sprintf('./Result/Result New/HV_GLHSSNAGO_Rand_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
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

function[hv,rt] = load_GLNAGO_DSS(ref,M,proDispName,selNum,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/DSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    rtDSS = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/GLHSSNAGO_DSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSelSet = load(fileName).dataSelSet;
    rt = load(fileName).runTime;
    rt = rtDSS+rt;
    fileName = sprintf('./Result/Result New/HV_GLHSSNAGO_DSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
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

function[hv,rt] = load_GLNAGO_IDSS(ref,M,proDispName,selNum,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/IDSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    rtIDSS = load(fileName).runTime;
    fileName = sprintf('./Result/Result New/GLHSSNAGO_IDSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    dataSelSet = load(fileName).dataSelSet;
    rt = load(fileName).runTime;
    rt = rtIDSS+rt;
    fileName = sprintf('./Result/Result New/HV_GLHSSNAGO_IDSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
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
    if exist(fileName,'file')
        hv = load(fileName).hv;
    else
        dataSelNorm = (dataSel-fmin)./(fmax-fmin);
        hv = stk_dominatedhv(dataSelNorm,ref);
        save(fileName,'hv');
    end
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