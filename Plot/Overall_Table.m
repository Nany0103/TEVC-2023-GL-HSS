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
    FID = fopen('./Table/HV_Overall_Comparison.tex','w');
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
        fprintf(FID,'\\hline \n');
    end
    fprintf(FID,'\\multicolumn{3}{c|}{(+/-)}');
    for colInd = 2:size(HV,2)
        fprintf(FID,'&(%d/%d)',sum(HVBet(:,colInd)),sum(HVWor(:,colInd)));
    end
    fprintf(FID,'\\\\ \n');
    fprintf(FID,'\\hline \n');
    fclose(FID);
    FID = fopen('./Table/RT_Overall_Comparison.tex','w');
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
                    fprintf(FID,'&%.2e (%s)',HVAve(rowInd,colInd),sym);
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
                    fprintf(FID,'&%.2e (%s)',HVAve(rowInd,colInd),sym);
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