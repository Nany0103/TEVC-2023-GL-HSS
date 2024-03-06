clc; clear;
obj = [3,5,8,10];
pro = {'Linear','Concave','Convex','I-Linear','I-Concave','I-Convex'};
numSol2 = 5000;
runNum = 21;
selNum = 100;
hvR = zeros(6,length(obj));
rtR = zeros(6,length(obj));

for proInd = 1:6
    proName = pro{proInd};
    hvAll = zeros(runNum,length(obj),2);
    rtAll = zeros(runNum,length(obj),2);
    for objInd = 1:length(obj)
        M = obj(objInd);
        H = getH(selNum,M);
        r = 1+1/H;
        ref = (1+1/H)*ones(1,M);
        for runInd = 1:runNum
            fprintf('%s_%d_%d\n',proName,selNum,runInd);
            % load pre time
            fileName = sprintf('./Result/Result New/TGAHSS_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat', ...
                M,proName,selNum,numSol2,runInd);
            rt_pre = load(fileName).runTime;
            %% load GL-HSS
            fileName = sprintf('./Result/Result New/GLHSSNAGO_TGA_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat', ...
                M,proName,selNum,numSol2,runInd);
            fileNameHV = sprintf('./Result/Result New/HV_GLHSSNAGO_TGA_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat', ...
                M,proName,selNum,numSol2,runInd);
            dataSelSet = load(fileName).dataSelSet;
            rt = rt_pre+load(fileName).runTime;
            if exist(fileNameHV,'file')
                hv = load(fileNameHV).hv;
            else
                hv = zeros(1,length(rt));
                for i=1:length(rt)
                    hv(i) = stk_dominatedhv(dataSelSet(:,:,i),ref);
                end
                parsave(fileNameHV,hv);
            end
            hvAll(runInd,objInd,1) = hv(end); rtAll(runInd,objInd,1) = rt(end);
            %% load LGI-HSS
            fileName = sprintf('./Result/Result New/LGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
                M,proName,selNum,runInd);
            fileNameHV = sprintf('./Result/Result New/HV_LGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
                M,proName,selNum,runInd);
            dataSel = load(fileName).dataSel;
            rt = load(fileName).runTime;
            if exist(fileNameHV,'file')
                hv = load(fileNameHV).hv;
            else
                hv = stk_dominatedhv(dataSel,ref);
                parsave(fileNameHV,hv);
            end
            hvAll(runInd,objInd,2) = hv; rtAll(runInd,objInd,2) = rt;
        end
    end
    %% calculate difference
    hvAve = zeros(2,length(obj));
    rtAve = zeros(2,length(obj));
    for i = 1:2
        hvAve(i,:) = mean(hvAll(:,:,i));
        rtAve(i,:) = mean(rtAll(:,:,i));
    end
    hvR(proInd,:) = hvAve(1,:)./hvAve(2,:);
    rtR(proInd,:) = rtAve(1,:)./rtAve(2,:);
end
%% plot hv figure
f = figure('Position',[200,200,750,450]);
hold on;
fs = 18;
lw = 1.5; ms = 15;
c = {'r','g','b','r','g','b'};
s = {'s','s','s','o','o','o'};
for proInd = [1,4,2,5,3,6]
    plot(1:4,hvR(proInd,:),['-',s{proInd}],'MarkerFaceColor',c{proInd},'MarkerSize',ms, ...
        'MarkerEdgeColor','k','LineWidth',lw,'Color',c{proInd});
end
ax = gca;
ax.XTick = 1:4;
ax.XTickLabel = {'3','5','8','10'};
ax.FontName = 'Times New Roman';
ax.FontSize = fs;
xlabel('Number of Objectives ({\itm})');
ylabel('Hypervolume Ratio');
ax.XGrid = 'on';
set(gca,'LineWidth', 2);
yTick = get(gca,'yTick');
yTickLabel = arrayfun(@(x) sprintf('%.3f',x), yTick,'uniformoutput',false);
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
saveas(f,'./Figure/Discussion_M/hv.emf');
close all;
%% plot runTime figure
f = figure('Position',[200,200,750,450]);
hold on;
c = {'r','g','b','r','g','b'};
s = {'s','s','s','o','o','o'};
for proInd = [1,4,2,5,3,6]
    plot(1:4,rtR(proInd,:),['-',s{proInd}],'MarkerFaceColor',c{proInd},'MarkerSize',ms, ...
        'MarkerEdgeColor','k','LineWidth',lw,'Color',c{proInd});
end
ax = gca;
ax.YScale = 'log';
ax.XTick = 1:4;
ax.XTickLabel = {'3','5','8','10'};
ax.YTick = [1e-3,1e-2,1e-1,1,10];
ax.YLim = [1e-3,10];
ax.FontName = 'Times New Roman';
ax.FontSize = fs;
xlabel('Number of Objectives ({\itm})');
ylabel('Runtime Ratio');
ax.XGrid = 'on';
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
saveas(f,'./Figure/Discussion_M/rt.emf');
close all;

function parsave(fileName,hv)
    save(fileName,'hv');
end

function H = getH(N,M)
    H = 1;
    while nchoosek(H+M,M-1) <= N
        H = H + 1;
    end
end