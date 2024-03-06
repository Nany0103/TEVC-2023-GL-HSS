% obj = [8,10,8,10,8,10,8,10,8,10,8,10,9,9,10];
% pro = {'DTLZ1','DTLZ1','DTLZ2','DTLZ2','MinusDTLZ1','MinusDTLZ1','MinusDTLZ2','MinusDTLZ2','WFG3','WFG3','DTLZ7','DTLZ7','RE91','DDMOP1','DDMOP4'};
% DPlus = [12,14,17,19,12,14,17,19,17,19,27,29,7,11,13];
% selNum = 100;
% runNum = 21;
% CandidateNum = zeros(1,length(pro));
% CandidateStd = zeros(1,length(pro));
% for proInd = 1:length(pro)
%     M = obj(proInd);
%     D = DPlus(proInd);
%     SA = zeros(1,runNum);
%     for runInd = 1:runNum
%         data = load(sprintf('./Data/EMOA/NSGAIII_%s_M%d_D%d_%d',pro{proInd},M,D,runInd)).UEA;
%         SA(runInd) = size(data,1);
%     end
%     CandidateNum(proInd) = mean(SA);
%     CandidateStd(proInd) = std(SA);
% end

write_table(CandidateNum,CandidateStd,pro,obj);

function write_table(Num,Std,pro,obj)
    FID = fopen('./Num_Of_Candidate_Set.tex','w');
    for rowInd = 1:6
        fprintf(FID,'%s',pro{(rowInd-1)*2+1});
        fprintf(FID,'&%.2e (%.2e)',Num((rowInd-1)*2+1),Std((rowInd-1)*2+1));
        fprintf(FID,'&%.2e (%.2e)',Num(rowInd*2),Std(rowInd*2));
        fprintf(FID,'\\\\ \n');
        fprintf(FID,'\\hline \n');
    end
    for rowInd = 13:length(Num)
        fprintf(FID,'%s',pro{rowInd});
        fprintf(FID,'&%.2e (%.2e)',Num(rowInd),Std(rowInd));
        fprintf(FID,'&');
        fprintf(FID,'\\\\ \n');
        fprintf(FID,'\\hline \n');
    end
    fclose(FID);
end