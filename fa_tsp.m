function fa_tsp()
%****************inputs*******************
global nFF;
global N;
global best;
global alpha;
global delta;

nFF = 50; %number of fireflies
movements = 30; %number of times a firefly moves
%gamma = 10; %light absorption coeffient
alpha = 0;
delta = 0;
iterations = 200; %number of times the FFs will evolve
file  = 'eil51.tsp'; %file name
minDist = 426;


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
    solutions(iteration) = best(1,N+1);        
end


[nodes1, nodes2] = getGraphNodes(best);

figure
plot (solutions);

figure
G = graph(nodes1,nodes2);
p = plot(G);
p.NodeColor = 'r';
p.XData = xValues;
p.YData = yValues;
% plot(G,'XData',xValues,'YData',xValues,'EdgeLabel');

disp('route');
disp(best(1:N));
disp('distance');
disp(best(1,N+1));
disp('difference from optimum solution');
disp(best(1,N+1) - minDist)




function setGlobalBest(val)
global best;
best = val;

function setGlobalGamma(val)
global gammaBest;
gammaBest = val;
    



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
        global N;
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
        FFs = zeros(FF,(nCity+2));
        for i= 1:FF
            FFs(i,1:nCity) = (randperm(nCity));
            FFs(i,(nCity+2)) = (rand()*10);
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
        global N;
        SortedFF = sortrows(FFs,(N+1));
        best = SortedFF(1,:);
        gamma = SortedFF(1,N+2);
        setGlobalBest(best);
        setGlobalGamma(gamma);


%********** Calculate Distance between 2 Fireflies (solutions) ****************

    function arcs = calDistSol(FF1, FF2)
        global N;
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
        global N;
        d  = FF(1,N+1);
        bri = 1/d;
    


%********** Calculate Attractiveness ****************

    function attr = calAttr(FF1, FF2)
        global N; 
        r = calDistSol(FF1, FF2);
        gamma = FF1(1, N+2);
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
        global N;
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

    function newFFs = newSols(FFset, movements, best) %TODO
        global N;
        nFF = size(FFset, 1);
        newFFs = zeros((nFF * movements)+1, N+2);
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
        
    function [nodes1, nodes2] = getGraphNodes(best)
        global N;
        nodes2 = zeros(N,1);
        for i = 1:(N-1)
            nodes2(i) = best(i+1);
        end
        nodes2(N) = best(1);
        nodes1 = best(1:N);

           
           
           
           
           
           
           
           
