function methods = getMethods(k,r,M)
    methods=[];
    iter = M;
    if k==1
        % Random partition (IR)
        for j=1:iter
            method = {};
            method.name = "random";
            method.seed = j*10+r;
            methods = [methods,method];
        end
    elseif k==2
        % Cosine similarity-based partition (IC)
        for j=1:iter
            method = {};
            method.name = "cosine";
            method.seed = j*10+r;
            methods = [methods,method];
        end
    elseif k==3
        % Objective value-based partition (IO)
         for j=1:iter
            method = {};
            method.name = "objective";
            method.seed = j*10+r;
            method.dim=mod(j-1,M)+1;
            methods = [methods,method];
         end   
    elseif k==5
        % Hybrid partition (IC + IO)
         for j=1:iter/2
            method = {};
            method.name = "cosine";
            method.seed = j*10+r;
            method.dim = 0;
            methods = [methods,method];
        end
         for j=1:iter/2
            method = {};
            method.name = "objective";
            method.seed = (iter/2+j)*10+r;
            method.dim = mod(j-1,M)+1;
            methods = [methods,method];
         end     
    elseif k==6
        % Hybrid partition (IO + IC)
         for j=1:iter/2
            method = {};
            method.name = "objective";
            method.seed = j*10+r;
            method.dim = mod(j-1,M)+1;
            methods = [methods,method];
         end   
         for j=1:iter/2
            method = {};
            method.name = "cosine";
            method.dim = 0;
            method.seed = (iter/2+j)*10+r;
            methods = [methods,method];
         end
    else
        error("No implementation!");
    end
end