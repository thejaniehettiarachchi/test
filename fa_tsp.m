function fa_tsp()
%****************inputs*******************

nFF = 25; %number of fireflies
movements = 20; %number of times a firefly moves
global gamma;
gamma = 0.01; %light absorption coeffient
iterations = 400; %number of times the FFs will evolve
file  = 'eil51.tsp'; %file name

%*************** Initialize variables ****************
xValues = 0; 
yValues = 0;
N = 0;
distMat = 0;
initFF = 0;
best = 0;
solutions = 0;



%**********     Read tsp file      **************

fileID = fopen(file);
%read file and store data in cell format
datacell = textscan(fileID, ' %f %f %f','CollectOutput',1);
fclose(fileID);

%convert data format from cell to matrix
datamat = cell2mat(datacell(1));

 
N = length(datamat); %No. of cities = No. of rows

cities = datamat(:,1); %city numbers
xValues = datamat(:,2); %x coordinates
yValues = datamat(:,3); %y coordinates


%****************** Start **************************


distMat = zeros(N,N);
solutions = zeros(iterations, 1);
distMat = disMat(distMat, xValues, yValues);
initFF = init(nFF, N);
newFF = calcObjFunc(initFF, N, distMat);
[newFF, best] = sort(newFF);

%**********     evolve      **************      
for iteration=1:iterations

    newPop = newSols(newFF, movements, best);
    newPop = calcObjFunc(newPop, N, distMat);
    newFF = selectFFs(newPop, nFF);
    disp(best.')
    solutions(iteration) = best(1,N+1)          
end
% figure
% plot (solutions);
disp (solutions);
% 
% disp('route');
% disp(best(1:N));
% disp('distance');
% disp(best(1,N+1));
% disp('difference from optimum solution');
% disp(best(1,N+1) - 426)

% function setGlobalG(val)
% global gamma
% gamma = val;
%     
% function G = getGlobalG
% global gamma
% r = gamma;


% N
% xValues
% yValues
% distMat
% initFF
%opt = [1 22 8 26 31 28 3 36 35 20 2 29 21 16 50 34 30 9 49 10 39 33 45 15 44 42 40 19 41 13 25 14 24 43 7 23 48 6 27 51 46 12 47 18 4 17 37 5 38 11 32];




%***************************************************
%                    FUNCTIONS
%***************************************************


%********** Create Distance Matrix ****************
    function distMat = disMat(distMat, xValues, yValues)
        N = size(xValues, 1);
        for i = 1:N
            for j = i:N
                %enter distance between i and j into distMat
                x1 = xValues(i); %x coodinate of c1
                y1 = yValues(i); %y coodinate of c1
                x2 = xValues(j); %x coodinate of c2
                y2 = yValues(j); %y coodinate of c2

                %calculate the distance between c1 and c2
                %dist = sqrt(((x1-x2)^2) + ((y1-y2)^2));
                X = [x1,y1;x2,y2];
                dist = pdist(X,'euclidean');
                distMat(i,j)= round(dist); %round off distance
                distMat(j,i) = distMat(i,j);
            end
        end

%********** Create Initial Population ****************
    function FFs = init(FF, nCity) %nFF = No. of FF, N = No. of cities
        FFs = zeros(FF,nCity);
        for i= 1:FF
            FFs(i,:) = (randperm(nCity));
        end
    


%********** Calculate Distance ****************
    function FFs = calcObjFunc(FFs, N, distMat) 
        num  = size(FFs, 1);
        for i = 1:num
            sum = 0;
            for j = 1:(N-1)
                x = (FFs(i,j));
                y = (FFs(i,(j+1)));
                d = distMat(x,y);
                sum = sum + d;
            end
            sum = sum + distMat(FFs(i, N),FFs(i, 1));
            FFs(i,(N+1)) = sum;
        end
%         disp(sum);
    



%********** Sort ****************

    function [SortedFF, best] = sort(FFs) 
        N = size(FFs.', 1)-1;
        SortedFF = sortrows(FFs,(N+1));
        FFs;
        SortedFF;
        best = SortedFF(1,:);
    


%********** Calculate Distance between 2 Fireflies (solutions) ****************

    function arcs = calDistSol(FF1, FF2)
         N = size(FF1.', 1)-1;
        arcs = 0;
        for i = 1:(N-1)
            temp1 = FF1(i);
            temp2 = FF1(i+1);
            j = 1;
            while FF2(j) ~= temp1
                j = j+1;
            end
            if j<N
                if FF2(j+1) ~= temp2
                    arcs = arcs + 1;
                end
            end
        end
        %arcs
    

%********** Calculate Brightness ****************

    function bri = brightness(FF)
        N = size(FF.', 1)-1;
        d  = FF(1,N+1);
        bri = 1/d;
    


%********** Calculate Attractiveness ****************

    function attr = calAttr(FF1, FF2)
        global gamma;
        r = calDistSol(FF1, FF2);
        attr0 = brightness(FF2);
        pow = -(gamma * r *r);
        attr = attr0 * exp(pow);
    

%********** Find More Attractive Firefly ****************

    function attrFF = getAttrFF (FF, FFset)
        nFF = size(FFset, 1);
        attrFF = -1;
        FFs = (randperm(nFF));
        for i= 1:nFF
            temp = FFset(FFs(i),:);
%             disp(brightness(FF));
%             disp(calAttr(FF, temp));
            if brightness(FF) < calAttr(FF, temp)
                attrFF = temp;
                break;
            end            
        end
    


%********** Inverse Mutation ****************

    function FF = invMutation(FF, length)
        N = size(FF.', 1)-1;
        startAt = randi(N, 1);
        temp = zeros(1,length);
        for i = 1:length
            if mod((startAt + (i-1)), N) ~= 0
                to  = mod((startAt + (i-1)), N);
            else
                to = N;
            end
            temp(i) = FF(to);
        end
        for i = 1:length
            if mod((startAt + (i-1)), N) ~= 0
                to  = mod((startAt + (i-1)), N);
            else
                to = N;
            end
            FF(to) = temp(length -(i-1));
        end
    


%********** Select new Solutions ****************

    function selectedFFs = selectFFs(FFs, nFF)
        FFs = sort(FFs);
        selectedFFs = FFs(1:nFF,:);
    


%********** New Solutions ****************

    function newFFs = newSols(FFset, movements, best)
        nFF = size(FFset, 1);
        N = size(FFset.', 1)-1;
        newFFs = zeros((nFF * movements)+1, N+1);
        for i = 1:nFF
            FF = FFset(i, :);
            attrFF = getAttrFF (FF, FFset);
            if attrFF == -1
                for j = 1:movements
                    move = randi([2,N],1);
                    newFFs(movements*(i-1) + j,:) = invMutation(FF, move);
                end
            else
                dist = calDistSol(FF, attrFF);
                if dist <2
                    dist = 2;
                end
                for j = 1:movements
                    move = randi([2,dist], 1);
                    newFFs(movements*(i-1) + j,:) = invMutation(FF, move);
                end
            end
        end
        newFFs((nFF * movements)+1, :) = best;
    
