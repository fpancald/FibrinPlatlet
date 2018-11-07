%this file makes a graph in the same way a hemGraph does, but in a
%cylindrical shape. No holes are placed in the network. This file is made
%to generate networks.

function Graph=RectangleGraphAlt() 
%Structural Origins of Fibrin Clot Rheology data
%fiberDensity = [50-100]
%branchDensity = [25-75]
clear all
close all
rng shuffle;
ErrorAngles=0.05;
ErrorLength=0.05;

fiberDensity = 1.0;
E1 = 100;%xlen 
E2 = 100;%ylen
E3 = 100;%zlen


totalFibers = E1 * E2 * E3 * fiberDensity;

maxDegree = 8;%unused, just use data

%biomaterials data
degree_nodes1 = 388;
degree_nodes2 = 0;
degree_nodes3 = 646;
degree_nodes4 = 128;
degree_nodes5 = 14;
degree_nodes6 = 0;
degree_nodes7 = 0;
degree_nodes8 = 0;
degree_nodes9 = 0;
degree_nodes10 = 0;
totalNodes = 1176;

degreeBiomat = [degree_nodes1,degree_nodes2,degree_nodes3,degree_nodes4,degree_nodes5,degree_nodes6,degree_nodes7,degree_nodes8]./totalNodes;
degreeChoice = degreeBiomat;


%% orientation of fibers data.
ErrorOrientation = 0.05; %percentage of edges allowed outside of 70% correct angle. 
prefOrientationAngle = 1.0; %gives average angle perference. 1.0 is z axis, 0.0 is x or y axis. 

orientation_bin_centers = [95,105,115,125,135,145,155,165,175,5,15,25,35,45,55,65,75,85];
no_flow_orientation = [0.052981146
0.052903428
0.053681837
0.052198755
0.058907136
0.061035482
0.05932573
0.056802632
0.054803359
0.058547209
0.058090158
0.059370446
0.062973519
0.063871322
0.051331607
0.050746768
0.046034342
0.04639513];

low_flow_orientation = [0.042720163
0.043731627
0.046191588
0.050547744
0.05092216
0.055036953
0.054090592
0.05904476
0.069506447
0.070157402
0.073671379
0.068156344
0.059806939
0.059408796
0.054537537
0.048074856
0.048625248
0.045649883];

high_flow_orientation = [0.031188673
0.032979896
0.029314074
0.027764853
0.028602296
0.034106976
0.047823461
0.072355511
0.133170597
0.143824643
0.109551209
0.062680509
0.054696462
0.049680018
0.038842595
0.038551434
0.034153912
0.03071293];

%sort data by key values
mapObjNoFlow = containers.Map(orientation_bin_centers,no_flow_orientation);
mapObjLowFlow = containers.Map(orientation_bin_centers,low_flow_orientation);
mapObjHighFlow = containers.Map(orientation_bin_centers,high_flow_orientation);

no_flow_orientationSort = cell2mat(values(mapObjNoFlow));
low_flow_orientationSort = cell2mat(values(mapObjLowFlow));
high_flow_orientationSort = cell2mat(values(mapObjHighFlow));
orientation_bin_centers_sort = cell2mat(keys(mapObjNoFlow));


%since the distribution should be symmetric, add data and reduce
orientation_bin_centers_sort_reduced = [];
no_flow_orientation_sort_reduced = [];
low_flow_orientation_sort_reduced = [];
high_flow_orientation_sort_reduced = [];
for i = 1:9
    orientation_bin_centers_sort_reduced(i) = orientation_bin_centers_sort(i);
    no_flow_orientation_sort_reduced(i) = no_flow_orientationSort(i) + no_flow_orientationSort(19-i);
    low_flow_orientation_sort_reduced(i) = low_flow_orientationSort(i) + low_flow_orientationSort(19-i);
    high_flow_orientation_sort_reduced(i) = high_flow_orientationSort(i) + high_flow_orientationSort(19-i);
end
%for now always use true. 
useOrientation = false;
AlignmentChoice='_';
orientationChoice = no_flow_orientation_sort_reduced;


if (orientationChoice == high_flow_orientation_sort_reduced)
    AlignmentChoice = 'high_flow_orientation';
end
if (orientationChoice == low_flow_orientation_sort_reduced)
    AlignmentChoice = 'low_flow_orientation';
end
if (orientationChoice == no_flow_orientation_sort_reduced)
    AlignmentChoice = 'no_flow_orientation';
end
%The choice is from suplementary materials
%The paper is :Structural basis for the nonlinear mechanics of fibrin networks under compression
%% Begin Parameters

%controls preferred edge length from paper
mu = .53;
sigma = .78; 

%for some reasaon, reversing the samples is working better. Normally it is
%giving too many large edges
sigmaGraphSample = 0.8;
muGraphSample = mu;% 0.8;
 
%prefDegree = 4; %minimal preferred
maxLength = 7;
densityEdge = 0.98; %don't set this to 1 or it will run forever

zLine = [0,0,1];


%% End Parameters
%num_edges = (degreeBiomat .* N .* [1,2,3,4,5,6,7,8])

N=round( totalFibers ./ sum(degreeBiomat.*[1,2,3,4,5,6,7,8]) );

%% generate random coordinates of nodes
%this makes a shape in quadrant 1 of length Len and radius Rad
%the center is (Len/2,Rad,Rad)
coord1=rand(1,N) * E1;
coord2=rand(1,N) * E2;
coord3 = rand(1,N) * E3;

%T1 holdes all node positions
T1=[cell(1,0),'nodes'];
for i=1:size(coord1,2)
    T1=[T1;[num2str(i),' x =   ',num2str(coord1(i)),'     y =    ',num2str(coord2(i)),'     z =    ',num2str(coord3(i))]];
end

%need to rewrite so no m
%instead of holding all nodes, it holds double the necesarry nodes. 
ro = N / (E1*E2*E3);
maxPossibleNbrs = ceil(8*(maxLength^3)*ro);%triple the average per square region
allLength = inf(length(coord1),maxPossibleNbrs);%holds lengths
allAngle = inf(length(coord1),maxPossibleNbrs);%angles close to 1 means z alignment, angles far from 1 mean no alignment
allId = zeros(length(coord1), maxPossibleNbrs);%holds id's

for i=1:length(coord1)
    for j=1:length(coord1)
        edge = [coord1(i)-coord1(j), coord2(i)-coord2(j),coord3(i)-coord3(j)];
        edgeLength = norm(edge);
        if (0<edgeLength) && (edgeLength<maxLength)
            %then we can place an edge. 
            for row = 1:maxPossibleNbrs
                if (allLength(i,row) == inf)
                    %then we have found a place to put our edge length
                    allLength(i,row) = edgeLength;
                    
                    %-1 is 90 degree, 1 is 0/180
                    allAngle(i,row) = 2*abs(dot(edge/edgeLength,zLine))-1; %unused right now
                    allId(i,row) = j;
                    break ;
                end
            end
        end
    end
end

%Now compute the density of nodes in a nbhd of radius maxLength around each point.
degreeNear = [];
for i = 1:length(coord1)
    num_nbrs_i = find(allLength(i,:) < maxLength);
    degreeNear(i) = length(num_nbrs_i);

end

%we sort the degree so the highest degree will be at the most dense place
[degreeAll,index] = sort(degreeNear);

%% Set the degree of each to the desired degree.
bins = floor(degreeChoice * N);

difference = N - sum(bins);
%choose three largest bins at random.
for i = 1:difference
    [~, ind] = sort(bins,'descend');
    choice = ind(ceil(3*rand));
    bins(choice) = bins(choice) + 1;
end

counter = 1;
for i = 1:length(degreeChoice)
    for j = 1:bins(i)
        %Set the degree
        degreeAll(counter) = i;%(prefDegree - 1) + i;
        counter = counter + 1;
    end
end

    
%% we construct the edges. Here choose the distribution for the edges. 


allPoints = [coord1;coord2;coord3]';

possibleEdges = [];
possiblePoints  = [];

T2 = [];

%choose random starting point

random = ceil(length(coord1)*rand);
startPoint = [coord1(random),coord2(random),coord3(random)];
indexStart = random;

%we will run the edge connecting until temp is empty
%we subtract 1 from the point where an edge is made 
%fill temp with the desired degree of each node computed above.
allEdges = zeros(length(coord1),1);
for i = 1:length(coord1)
    allEdges(index(i)) = degreeAll(index(i));
end
sum_temp = sum(allEdges);


%% Making The Graph 
while (densityEdge * sum_temp/2 > length(T2))
    %each end point is weighted
    possiblePoints = nonzeros(allId(indexStart,:));
    probOfStep = [];
    
    for i = 1:length(possiblePoints)   
        length_i = [allLength(indexStart, i)];%lengths associated to id
        %angle_i = 0.5*[allAngle(indexStart, i)];
        %%%%%
        probOfStep(i) = lognpdf(length_i,muGraphSample,sigmaGraphSample);%weighted
        
        %re-wieght using angles. 
        if (useOrientation)
            %probOfStep(i) = probOfStep(i) * angle_i;
        end
        
    end
       
    %choose the endpoint according to this weight.
    indexChoice = find(probOfStep == max(probOfStep));
    prefEdgeIndex = possiblePoints(indexChoice);
    %[row, col] = find(allId(indexStart,:) == prefEdgeIndex)
    prefEdgeLen = allLength(indexStart, indexChoice);
 
    %run this part until a endpoint is chosen each iteration attempts to
    %place a single edge.
    
    counter = length(T2)+1;
    while length(T2) < counter
        
        %Initial choice test. the simple line worked better
        %width = max(probOfStep) - min(probOfStep);
        %random = rand  - (1-prefEdgeLen) + width;
        %draw = lognpdf(random, mu, sigma);
        draw = lognrnd(muGraphSample,sigmaGraphSample);
        
        minProb =(min(abs(probOfStep-draw)));
        %probOfStep

        choice = find(abs(probOfStep-draw) == minProb,1);
        %choice
        indexEnd = possiblePoints(choice);
        endPoint = allPoints(indexEnd,:);
        
%% Duplicates: If there is already an edge, lia == (true, true). 
%while this holds, we check all possible edges and then break when there
%are no more. 
        lia = [0,0];
        for i = 1:size(T2,1)
            if ((indexStart == T2(i,1)) && (indexEnd ==  T2(i,2)))...
                        || ((indexStart == T2(i,2)) && (indexEnd == T2(i,1)))
                lia(1,1) = true;
                lia(1,2) = true;
            end
        end

        %to avoid duplicates 
        while (lia(1,1) == true) && (lia(1,2) == true)%check if the edge is already in T2
            %re-choose the end point
            probOfStep(choice) = 0;
            %random = rand  - (1-prefEdgeLen) + width;
            %draw = lognpdf(random, mu, sigma); %+ min(probOfStep);
            draw = lognrnd(muGraphSample,sigmaGraphSample);
            minProb =(min(abs(probOfStep-draw)));
            choice = find((abs(probOfStep-draw) == minProb),1);
            indexEnd = possiblePoints(choice);

            endPoint = allPoints(indexEnd,:);
            
            %make sure you reset the lia
            lia = [0,0];
            for i = 1:length(T2)
                if (((indexStart == T2(i,1)) && (indexEnd ==  T2(i,2)))...
                        || ((indexStart == T2(i,2)) && (indexEnd == T2(i,1))))
                    lia(1,1) = true;
                    lia(1,2) = true;
                end
            end
            
            if any(probOfStep ~= 0)
                continue;
            else
                %if you reach here, all probOfStep will be zero.
                break
            end
        end
        
%% Adding Edges is done on a first needed basis. If allEdges has more room and probOfstep is not all zeros 
% a new edge is created.
         %indexEnd
         %endPoint
        if (allEdges(indexEnd) > 0) &&...
                (any(probOfStep ~= 0) &&...
                (true ))
            T2 = [T2; [indexStart, indexEnd]];
            allEdges(indexStart) = allEdges(indexStart) - 1;%account for created edge by subtracting off
            allEdges(indexEnd) = allEdges(indexEnd) - 1;
            
            %indexStart = indexEnd; %use this for random walk constructor
            %startPoint = endPoint; %use this for random walk
            
            %next point is chosen by mean
            possibleStarts = []; %holds top 
            threshold = mean(allEdges); %notice this is different from the bottom
            for i = 1:length(allEdges)
                if (allEdges(i)> threshold)
                    possibleStarts = [possibleStarts; i];
                end
            end
            
            random = ceil(length(possibleStarts)*rand);
            newIndex = possibleStarts(random);
            startPoint = [coord1(newIndex),coord2(newIndex),coord3(newIndex)];
            indexStart = newIndex;
            
            %now that we have a new index, we need to regenerate the
            %probabilities of possible end points.
            possiblePoints = nonzeros(allId(indexStart,:));
            probOfStep = [];

            for i = 1:length(possiblePoints)   
                length_i = [allLength(indexStart, i)];
                %angle_i = 0.5*[allAngle(indexStart, i) ];
                probOfStep(i) = lognpdf(length_i,muGraphSample,sigmaGraphSample);%weighted
                
                %re-wieght using angles. 
                if (useOrientation)
                %    probOfStep(i) = probOfStep(i) * angle_i;
                end
            end
        %if
        elseif any(probOfStep ~= 0)
            %if the choice was not chosen, we set its prob to zero
            %so it will not be chosen again.
            probOfStep(choice) = 0;
        else
            %you only reach this point if you cannot form an edge and all
            %probabilities of endpoints are zero.
            %chose a new point from top ~50% of allEdges
            
            possibleStarts = []; %holds top ~50% 
            threshold = mean(allEdges);
            for i = 1:length(allEdges)
                if (allEdges(i)> threshold)
                    possibleStarts = [possibleStarts; i];
                end
            end
            
            random = ceil(length(possibleStarts)*rand);
            newIndex = possibleStarts(random);
            startPoint = [coord1(newIndex),coord2(newIndex),coord3(newIndex)];
            indexStart = newIndex;
            %force the next choice to have an available edge in allEdges
            %vector
            while (allEdges(indexStart) == 0)
                random = ceil(length(possibleStarts)*rand);
                newIndex = possibleStarts(random);
                startPoint = [coord1(newIndex),coord2(newIndex),coord3(newIndex)];
                indexStart = newIndex;
            end
            break
        end
        
    end
    counter = counter+1; %this makes the while loop run until an end is chosen
end 
    

%calculate degrees
degrees = zeros(N,1);
for i = 1:N
    degrees(i) = sum((sum(i == T2)));
end



%% Perterbation of points to fit edge lengths

lengths = [];
currentAngles = [];
for i=1:size(T2,1)
    x1 = coord1(T2(i,1));
    x2 = coord1(T2(i,2));
    y1 = coord2(T2(i,1));
    y2 = coord2(T2(i,2));
    z1 = coord3(T2(i,1));
    z2 = coord3(T2(i,2));
    edgeTemp = [x1 - x2, y1-y2, z1-z2];
    normTemp = norm(edgeTemp);
    lengths(i) = normTemp; 
    currentAngleTemp = 2*abs(dot(zLine, edgeTemp/normTemp))-1;
    currentAngles(i) = (currentAngleTemp*(-45) + 45);
   
end 
binSize = 0.1;
numSamples = length(lengths);%number of edges
lognSamples = zeros(numSamples,1);
tempMaxLength = maxLength+1;
for i =1:numSamples
    sample = tempMaxLength;
    while (sample >= tempMaxLength)
        sample = lognrnd(mu,sigma);
        lognSamples(i) = sample;
    end
end
BMIN = 0;
BMAX = maxLength+2.0; %safety feature
boundMin = min(min(coord1),0.0);
boundMax = max(coord1);
hist = histogram(lognSamples,'binwidth',binSize,'BinLimits',[BMIN,BMAX]);
preferredNumEdgesPerBin = hist.Values;
binEdges = hist.BinEdges;
[parmhat] = lognfit(lognSamples, 'lognormal'); %returns parameters for lognormal dist fitting data.


hist1 = histogram(lengths,'binwidth',binSize,'BinLimits',[BMIN,BMAX]);
currentNumEdgesPerBin = hist1.Values;

beginErrorEdges = length(T2);%sum(abs(preferredNumEdges - currentNumEdges));
errorCurrentEdges = beginErrorEdges;

%% we are going to scale the degrees  based on the left and right side
preferredNumAnglesPerBin = orientationChoice*length(T2);

angleBinDist = 5;
anglebinEdges = [orientation_bin_centers_sort_reduced-angleBinDist,90];
TempHist = histogram(discretize(currentAngles,anglebinEdges));
currentNumAnglesPerBin=TempHist.Values;

%we'll updatte lenthsUpdate and anglesUpdate
lengthsUpdate = lengths';
anglesUpdate = currentAngles';
maxAngleError = ErrorAngles*beginErrorEdges;
%
while ((errorCurrentEdges> ErrorLength * beginErrorEdges) || (errorCurrentAngles> ErrorAngles*beginErrorEdges))
  %  && (errorCurrentAngles > ErrorOrientation * beginErrorAngles))
    pointChoice = 1;
    testrnd = rand;
    
    %80% choose maximum error bins, 20% choose random
    if (rand < 0.6)
        maxError = max((currentNumEdgesPerBin - preferredNumEdgesPerBin));
        mostErrorBins = find((currentNumEdgesPerBin - preferredNumEdgesPerBin) > 0.7 * maxError);
        %the bins in mostErrorBins should be given the most attention we'll use
        %lengthsUpdate to choose the next pointChoice. 
        msize = numel(mostErrorBins);
        binChoice = mostErrorBins(randperm(msize, 1));
        edgeLenMin = (binChoice-1)*binSize;
        edgeLenMax = (binChoice)*binSize;
        [a,b] = find((lengthsUpdate>edgeLenMin) & (lengthsUpdate<edgeLenMax));
        row = randperm(length(a),1);
        col = randperm(2,1);
        pointChoice = T2(row,col);
    else
        pointChoice = ceil(rand*length(coord1));
    end

    
    %instead of choosing a random point, choose a point from top 20% error
    %bars.
    edgesAreInRange = 0;
    while (edgesAreInRange == 0)
        %choose perturbation 
        nodeChangeIsInRange = 0;
        while (nodeChangeIsInRange == 0)
            %choose node in range
             
            changeX = E1/4*(rand - 0.5);
            changeY = E2/4*(rand - 0.5);
            changeZ = E3/4*(rand - 0.5);
            
            prefOrientationAngle=0.71;%tester
            if (useOrientation)
                if (prefOrientationAngle > 0.7)
                    changeZ = E3 *(rand - 0.5);
                elseif ((prefOrientationAngle < 0.7) && (prefOrientationAngle > 0.5))
                    changeZ = 0.5*E3 *(rand - 0.5);
                else
                    changeZ =0.25*E3 *(rand - 0.5);
                end
            end

            xTestPos = coord1(pointChoice) + changeX;
            yTestPos = coord2(pointChoice) + changeY;
            zTestPos = coord3(pointChoice) + changeZ;
            distMax = max(max(xTestPos,yTestPos),zTestPos);
            distMin = min(min(xTestPos,yTestPos),zTestPos);
            if ((distMax < boundMax) && (distMin>boundMin))
                nodeChangeIsInRange = 1;
            end
        end

        [x,y] = (find(T2 == pointChoice));
        
        %get nbrs
        neighbors = zeros(length(x),1);
        if (length(neighbors) == 0)
            pointChoice = ceil(rand*length(coord1));
        end
        for nbr = 1:length(x)
            row = x(nbr);
            col = y(nbr);
            nbrCol = 0;
            if (col == 2)
                nbrCol = 1;
            else
                nbrCol = 2;
            end
            neighbors(nbr) = T2(row,nbrCol);
        end
        
        degree = length(neighbors);
        lengthCurrent = zeros(degree,1);
        lengthTest = zeros(degree,1);
        angleCurrent = zeros(degree,1);
        angleTest = zeros(degree,1);
        
        for edge = 1:degree
            %calc edges
            nbr = neighbors(edge);
            currentEdge = [coord1(pointChoice)-coord1(nbr),...
                coord2(pointChoice)-coord2(nbr),...
                coord3(pointChoice)-coord3(nbr)];
            
            lengthCurrent(edge) = norm(currentEdge);
            
            testEdge = [xTestPos-coord1(nbr),...
                yTestPos - coord2(nbr),...
                zTestPos - coord3(nbr)];
            
            lengthTest(edge) = norm(testEdge);

            %0-90
            tempAngleCurrent = 2 * abs(dot(zLine,currentEdge/norm(currentEdge) )) - 1;
            angleCurrent(edge) = (tempAngleCurrent*(-45) + 45);
            
            tempAngleTest = 2 * abs(dot(zLine,testEdge/norm(testEdge) )) - 1;
            angleTest(edge) = (tempAngleTest*(-45) + 45);
            
        end
        
        if (max(lengthTest) < BMAX)
            edgesAreInRange = 1;
        end
    end
    

    
    
    binCurrentLength = discretize(lengthCurrent,binEdges);
    binTestLength = discretize(lengthTest,binEdges);
    binCurrentAngle = discretize(angleCurrent, anglebinEdges);
    binTestAngle = discretize(angleTest, anglebinEdges);
    
    %impose on test histogram.
    %now we have test lengths vs original lengths. Compare them and see if the
    %histogram is more accurate. If so, accept the movement. Else continue. 

    testCurrentNumEdgesPerBin = currentNumEdgesPerBin;
    
    errorCurrentEdges = sum(abs(currentNumEdgesPerBin - preferredNumEdgesPerBin));
    for num = 1:degree
        binCurrent = binCurrentLength(num);
        binTest = binTestLength(num);
        if isnan(binTest)
            disp('binTest out of bounds');
        end
        testCurrentNumEdgesPerBin(binTest) = testCurrentNumEdgesPerBin(binTest)+1;
        testCurrentNumEdgesPerBin(binCurrent) = testCurrentNumEdgesPerBin(binCurrent) - 1;
    end    
    errorTestLengths = sum(abs(testCurrentNumEdgesPerBin - preferredNumEdgesPerBin));

    %now do the same for angles
   % TempHist1 = histogram(discretize(anglesUpdate,anglebinEdges));
   % currentNumAnglesPerBin=TempHist1.Values;

    testCurrentNumAnglesPerBin = currentNumAnglesPerBin;
    if (useOrientation)
        errorCurrentAngles = sum(abs(currentNumAnglesPerBin - preferredNumAnglesPerBin));
    else 
        errorCurrentAngles=0;
    end
    
    for num = 1:degree
        binCurrentTempAngle = binCurrentAngle(num);
        binTestTempAngle = binTestAngle(num);
        if isnan(binTestTempAngle)
            disp('binTest out of bounds');
        end
        testCurrentNumAnglesPerBin(binTestTempAngle) = testCurrentNumAnglesPerBin(binTestTempAngle)+1;
        testCurrentNumAnglesPerBin(binCurrentTempAngle) = testCurrentNumAnglesPerBin(binCurrentTempAngle) - 1;
        if (sum(testCurrentNumAnglesPerBin<0)>0)
            
        end
    end    
    errorTestAngles = sum(abs(testCurrentNumAnglesPerBin - preferredNumAnglesPerBin));

    
    %final perturbation acceptence
    condition=(errorTestLengths <= errorCurrentEdges);
    conditionOrientation=((errorTestLengths <= errorCurrentEdges) && ((errorTestAngles <= errorCurrentAngles) || (errorCurrentAngles<maxAngleError)) );
    if (useOrientation)
        condition=conditionOrientation;
    end
    
    if (condition)        %accept change and update point and histogram
        %update angle count
       % changeInAngles = errorTestAngleTemp - errorCurrentAngleTemp;
       % errorCurrentAngles = errorCurrentAngles - changeInAngles;
        
        coord1(pointChoice) = xTestPos;
        coord2(pointChoice) = yTestPos;
        coord3(pointChoice) = zTestPos;
        currentNumEdgesPerBin = testCurrentNumEdgesPerBin;    
        currentNumAnglesPerBin = testCurrentNumAnglesPerBin;
        %update lengthsUpdate and anglesUpdate to help choose new points
        %update using lengthTest and angleTest
        [a,b] = find(T2 == pointChoice);
        for i = 1:size(a,1)
            lengthsUpdate(a(i)) = lengthTest(i);
            anglesUpdate(a(i)) = angleTest(i); 
        end
    end
    
end


Fixed = 0;
all=[coord1',coord2',coord3'];
Graph.nd=all;
Graph.ed=T2;
Graph.Data={Fixed};

DateString=datestr(clock);
DateString(DateString==':')='_';
DateString(end-8:end)=[];

DateString=['Rec_',num2str(ro),'_', num2str(E2), '_', num2str(E2), '_', num2str(E3),'_',num2str(7)];
if (useOrientation == true)   
    DateString=['Rec_',num2str(ro), '_', num2str(E2), '_', num2str(E2), '_', num2str(E3),'_',AlignmentChoice];
end
Graph.String=DateString;
savdir = '/afs/crc.nd.edu/user/s/sbritto2/networkgen/MatlabFiles';
 save(fullfile(savdir,[DateString,'.mat']),'Graph');
 