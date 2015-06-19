function OrderID

    clear;
    TRUE = uint8(1);
    FALSE = uint8(0);

    numberOfSwaps = uint16(0); % The maximum number of pairs to swap in the section
    maxIterations = 5000; % How many synthetic shuffled sections to calculate to define the random distribution
    dataFacies = zeros(1,1);
    dataThick = zeros(1,1);
    markovOrderMetric = zeros(1,1);
    runsOrderMetric = zeros(1,1);
    maxRun = 3;
    minRun = 0;
    runBinIncrement = 0.05;
    runRange = maxRun - minRun;
    shuffledFacies = zeros(1,1);
    shuffledThick = zeros(1,1);
    dataShuffledMultiMarkovOrderMetric = zeros(1,maxIterations);
    multiRunsOrderMetric = zeros(1,maxIterations);

    % Read in the section facies and thickness data - you might need to edit the directory names
    fileName = input('Enter succession data path and file name (without any .dat or .txt extension)\ne.g. ../outcropSections/test100 ../syntheticSections/syntheticOrdered100units10facies\nFilename is : ','s');
    fullFileName = strcat(fileName, '.txt');

    if exist(fullFileName, 'file')
        % Read section data, thickness and facies
        dataIn = load(fullFileName);
        dataThick = dataIn(:,1);
        dataFacies = dataIn(:,2);
    else
        fprintf('Succession data file %s does not exist. Please run this script again and enter a valid path and filename\n', fullFileName);
        return;
    end

    % Read colour map file, one colour and one label text per facies
    fullColourMapFileName = strcat(fileName, 'LithoCol.txt');
    if exist(fullFileName, 'file')
        fid = fopen(fullColourMapFileName);
        while ~feof(fid)
            file = textscan(fid, '%u8 %f %f %f %s');
        end
    else
        fprintf('Colour code and labels file %s does not exist. Please run this script again and enter a valid path and filename\n', fullFileName);
        return;
    end

    ticks = file {5}();
    for i=1:4
        faciesColours(:,i)= double(file{i}());
    end
    fclose(fid);

    % Calculate the basic stats on the dataFacies array
    sectionLength = max(size(dataFacies));
    maxNumbOfFacies = max(dataFacies);
    sectionMean = mean(dataThick);
    fprintf('For %d units, total %d facies, mean unit thickness %4.3f m\n', sectionLength, maxNumbOfFacies, sectionMean);
    % Now use these data to define the number of swaps needed and dimension results array accordingly
    numberOfSwaps = sectionLength;

    % Calculate and output the order metric for the entire data succession
    markovOrderMetric = calculateTPMatrixAndOrderMetric(dataFacies, 0);  % note the zero is the TP matrix plot flag, set to false to not trigger a plot
    runsOrderMetric = calculateRunsOrderMetric(dataThick);
    fprintf('Markov metric for strata is %4.3f\nRuns analysis metric for strata is %5.4f\n', markovOrderMetric, runsOrderMetric);

    % Now calculate the metrics for many iterations of a random model

    % Shuffle the observed section and calculate the order metric and DW stat each time
    for j = 1:maxIterations;
        [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(dataFacies, dataThick, numberOfSwaps);
        multiMarkovOrderMetricDataShuffled(j) = calculateTPMatrixAndOrderMetric(shuffledFacies, 0); % note the zero is the TP matrix plot flag, set to false so no plot
        multiRunsOrderMetricDataShuffled(j) = calculateRunsOrderMetric(shuffledThick); 
    end

    % Stats on the shuffled section random model
    meanMultiMarkovDataShuffled = mean(multiMarkovOrderMetricDataShuffled);
    stdDevMultiMarkovDataShuffled = std(multiMarkovOrderMetricDataShuffled);
    bins = 0:0.02:1.00;
    multiMarkovOrderMetricDataShuffledHistData=histc(multiMarkovOrderMetricDataShuffled, bins) / maxIterations; % Calculate frequency bins with histc but / by iterations to give relative freq
    markovPValueSum = sum(multiMarkovOrderMetricDataShuffledHistData(int16(markovOrderMetric*50:length(bins)))); % area under curve from m to max m value 1

    meanMultiRunsDataShuffled = mean(multiRunsOrderMetricDataShuffled);
    stdDevMultiRunsDataShuffled = std(multiRunsOrderMetricDataShuffled);
    bins = minRun:runBinIncrement:maxRun;
    multiRunsOrderMetricDataShuffledHistData = histc(multiRunsOrderMetricDataShuffled, bins)/maxIterations; % Calculate frequency bins with histc but / by iterations to give relative freq
    runBinIndex = 1 + int16(((runsOrderMetric-minRun)/runRange)*(runRange/runBinIncrement)); % Position of runs stat in the histogram
    runsPValueSum = sum(multiRunsOrderMetricDataShuffledHistData(runBinIndex:length(multiRunsOrderMetricDataShuffledHistData))); % area under curve from r to max run value

    fprintf('For %d iterations of a SHUFFLED DATA model\n', maxIterations);
    fprintf('Markov stats mean %5.4f std dev %5.4f Markov order metric P Value %5.4f\n', meanMultiMarkovDataShuffled, stdDevMultiMarkovDataShuffled, markovPValueSum);
    fprintf('Runs stats mean %5.4f std dev %5.4f Runs analysis metric P Value %5.4f\n', meanMultiRunsDataShuffled, stdDevMultiRunsDataShuffled, runsPValueSum);

    % Plot the results
    scrsz = get(0,'ScreenSize'); % screen dimensions vector

    % Plot the original and shuffled vertical sections, alongside the facies and thickness
    % order plots showing the outcrop datapoint and the multiple iteration shuffled section
    % metric frequency distribution
    f2 = figure('Visible','on','Position',[1 scrsz(4)/4 (scrsz(3)/3)* 2.0 (scrsz(4)/3)*2]);

    % Subplot for the data vertical section, plotted on the left
    subplot('Position',[0.04 0.1 0.05 0.85]);
    hold on
    cumThick = 0;
    for j=1:sectionLength
        fCode = dataFacies(j);
        yco = [cumThick, cumThick, cumThick + dataThick(j), cumThick + dataThick(j)];
        xco = [0, dataFacies(j), dataFacies(j), 0];
        faciesCol = [faciesColours(fCode,2) faciesColours(fCode,3) faciesColours(fCode,4)];
        patch(xco, yco, faciesCol,'EdgeColor','none');
        cumThick = cumThick + dataThick(j);
    end

    grid on;
    set(gca,'Layer','top');
    xlim([0 maxNumbOfFacies])
    xlabel('Facies code');
    ylabel('Thickness (m)');
    set(gca,'XTick', 1:(maxNumbOfFacies),'XTickLabel', ticks, 'TickDir', 'out');

    % Subplot for the Markov order analysis histogram   
    subplot('Position',[0.42 0.57 0.55 0.38]);
    hold on
    bins = 0:0.02:1.00; % Make sure bins is set correctly for Markov plots
    bar(bins, multiMarkovOrderMetricDataShuffledHistData, 'EdgeColor','none', 'BarWidth', 1, 'FaceColor',[0.2 0.4 0.7]); % Colour is dark slate blue
    maxFreq = max(multiMarkovOrderMetricDataShuffledHistData) *1.1; % This is needed to scale the plot
    x = [markovOrderMetric markovOrderMetric];
    y = [0 maxFreq];
    line(x,y, 'color', [0.80 0.00 0.00], 'linewidth', 3.0); % Colour is dark red
    grid on;
    axis([0 1 0 Inf]);
    set(gca,'Layer','top');
    xlabel('Markov Order Metric for Facies');
    ylabel('Relative Frequency');

    % Subplot for the runs analysis histogram    
    subplot('Position',[0.42 0.1 0.55 0.38]);
    hold on
    bins = minRun:runBinIncrement:maxRun; % Make sure that bins is set correctly for Runs analysis plots
    bar(bins, multiRunsOrderMetricDataShuffledHistData, 'EdgeColor','none', 'BarWidth', 1, 'FaceColor',[0.2 0.4 0.7]);
    maxFreq = max(multiRunsOrderMetricDataShuffledHistData) * 1.1; % This is needed to scale the plot

    lineColor = [0.80 0.00 0.00];
    x = [runsOrderMetric runsOrderMetric];
    y = [0 maxFreq]; % Draw the data line from y=0 to y=max frequency of the three histograms
    line(x,y, 'color', lineColor, 'linewidth', 3.0);

    grid on;
    axis([0 Inf 0 Inf]);
    set(gca,'Layer','top');
    xlabel('Runs Analysis Order Metric for Thickness ');
    ylabel('Relative Frequency');

    % Plot the histogram of the facies frequencies
    subplot('Position',[0.13 0.1 0.23 0.38]);
    hold on
    frequencyDataFacies=histc(dataFacies,1:maxNumbOfFacies);
    for j=1:maxNumbOfFacies
        xco = [j-1, j-1, j, j];
        yco = [0, frequencyDataFacies(j), frequencyDataFacies(j), 0];
        faciesCol = [faciesColours(j,2) faciesColours(j,3) faciesColours(j,4)];
        patch(xco, yco, faciesCol,'EdgeColor','black');
    end
    set(gca,'XTick', 0.5:(maxNumbOfFacies-0.5),'XTickLabel', ticks);
    xlabel('Lithofacies');
    ylabel('Frequency');

    % plot the TP matrix
    subplot('Position',[0.13 0.57 0.23 0.38]);
    hold on
    calculateTPMatrixAndOrderMetric(dataFacies, 1);
    axis tight;

end

function markovOrderMetric = calculateTPMatrixAndOrderMetric(faciesSect, plotFlag)

    % find the number of elements in the succession. NB size returns both
    % matrix dimensions so max ensures nz is the biggest which should be
    % the length of the section
    nz = max(size(faciesSect));
    
    % Find the maximum facies code used in the facies succession - this is the
    % size for both dimensions of the TP matrix which can now be defined
    nFacies = max(faciesSect);
    m = 0; % number of transitions, not really needed but never mind
    TFMatrix = zeros(nFacies, nFacies);
    TPMatrix = zeros(nFacies, nFacies);

    % Now loop through the elements in the succession and for each different facies from-to transition,
    % increment the appropriate cell in the matrix
    for i =1 : nz-1
        fromFacies = faciesSect(i);
        toFacies = faciesSect(i+1);
        % mark transitions between different facies
        if fromFacies > 0 && toFacies > 0 && fromFacies ~= toFacies % Make sure facies codes are not zero because zero values would record an error
            TFMatrix(fromFacies, toFacies) = TFMatrix(fromFacies, toFacies) + 1; % increment the appropriate value in the tp matrix
            m = m + 1;
        end     
    end

    %Now calculate the transition probability matrix from the transition frequency
    %matrix
    rowSums=sum(TFMatrix,2); % Calculates the sum of each row in TF matrix and stores as vector rowSums
    for k=1:nFacies
        for j=1:nFacies
            if rowSums(k) > 0 % if rowsum > 0 divide TF value by row sum to get transition probability
                TPMatrix(k,j)=TFMatrix(k,j) / rowSums(k);
            else
                TPMatrix(k,j) = 0;
            end
        end
    end
    
    % Now calculate the Markov order metrics for the maximum diagonal
    diagMetric = zeros(1,nFacies-1);
    for j=1:nFacies-1;
        diagMetric(j) = (sum(diag(TPMatrix,j)) + sum(diag(TPMatrix,-(nFacies-j))) )/nFacies; % calculate a mean for each FULL diagonal in the matrix
    end
    markovOrderMetricDiags = max(diagMetric)- min(diagMetric);
    markovOrderMetric = markovOrderMetricDiags;
    
    if plotFlag == 1
        
        TPMatrix
        a1 = gca;
        TPCellSize = 1.0;
        matrixTopYco = (nFacies + 1) * TPCellSize;
        matrixBottYco = 0.0;
        
        for i = 1:nFacies
            for j = 1:nFacies
                yco = [matrixBottYco+((j-0.5)*TPCellSize) matrixBottYco+((j-0.5)*TPCellSize) matrixBottYco+((j+0.5)*TPCellSize) matrixBottYco+((j+0.5)*TPCellSize)];
                xco = [(i-0.5)*TPCellSize (i+0.5)*TPCellSize (i+0.5)*TPCellSize (i-0.5)*TPCellSize];
                labelStr = sprintf('%3.2f',TPMatrix(j,i));
                
                if TPMatrix(j,i) <= 0.5
                    colVect = [1 TPMatrix(j,i)/0.5 TPMatrix(j,i)/1.1111111]; % gradation of colours from red (p=0) to orange-yellow (p=0.5)
                else
                    colVect = [1-(TPMatrix(j,i)) 1 1-(TPMatrix(j,i)/1.111111)]; % gradation of colours from orange-yellow (p=0.5) to green (p=1)
                end
               
                if i==j
                    colVect = [0.9 0.9 0.9]; % Matrix diagonal is a special case so colour light grey
                end
                
                patch(xco, yco, colVect);
                
                if i ~= j % So do not add text labels to the diagonal
                  text(double(((i-0.7)*TPCellSize))+(TPCellSize*0.33), double(matrixBottYco+(j*TPCellSize)-0.3), labelStr,'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
                end
            end
        end
        
        set(a1,'XTick',1:nFacies);
        set(a1,'YTick',1:nFacies);
        xlabel('Facies Code: To', 'FontSize',10);
        ylabel('Facies code: From', 'FontSize',10);
        
    end
end

function runsOrderMetric = calculateRunsOrderMetric(thicknesses)

    % find the number of units in the succession and declare arrays accordingly
    nz = max(size(thicknesses));
    deltaThick = zeros(1,nz);
    runsUp = zeros(1,nz);
    runsDown = zeros(1,nz);

    % Calculate the change in thickness between successive units
    i = 1:nz-1;
    j =2:nz; % so j = i + 1 therefore thickness change is thickness(j) - thickness(i)
    deltaThick(i) = thicknesses(j) - thicknesses(i);

    if deltaThick(1) > 0 runsUp(1) = 1; end
    if deltaThick(1) < 0 runsDown(1) = 1; end

    for i=2:nz
        if deltaThick(i) > 0 runsUp(i) = runsUp(i-1)+1; end
        if deltaThick(i) < 0 runsDown(i) = runsDown(i-1)+1; end
    end

    runsUpNormSum = (sum(runsUp)/nz);
    runsDownNormSum = (sum(runsDown)/nz);
    runsOrderMetric = (runsUpNormSum + runsDownNormSum);
end

function [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(sectFacies, sectThick, totalSwaps)
% function to shuffle the facies succession to ensure a random configuration

    % Make copies of the original data in new arrays that will be used to
    % store the shuffled sections
    shuffledFacies = sectFacies;
    shuffledThick = sectThick;
    n = uint16(max(size(shuffledFacies)));
    j = 0;
    infiniteStopCount = 0;
    
    while j < totalSwaps && infiniteStopCount < 1000000
        
        % Select two unit numbers randomly to be swapped
        unit1 = uint16((rand * (n-1)) + 1);
        unit2 = uint16((rand * (n-1)) + 1);
        
        % Need to check above and below for both positions that swapping will not put same
        % facies adjacent to one another and cause a transition to self
        swapFacies1 = shuffledFacies(unit1);
        if unit1 > 1 swapFacies1Below = shuffledFacies(unit1-1); else swapFacies1Below = 0;end
        if unit1 < n swapFacies1Above = shuffledFacies(unit1+1); else swapFacies1Above = 0;end
        
        swapFacies2 = shuffledFacies(unit2);
        if unit2 > 1 swapFacies2Below = shuffledFacies(unit2-1); else swapFacies2Below = 0;end
        if unit2 < n swapFacies2Above = shuffledFacies(unit2+1); else swapFacies2Above = 0;end
        
        % So compare facies in their new positions with the facies above and below and
        % only swap and increment loop counter if NOT the same...
        if swapFacies1Below ~= swapFacies2 && swapFacies1Above ~= swapFacies2 && swapFacies2Below ~= swapFacies1 && swapFacies2Above ~= swapFacies1
            
            %Swap the facies
            temp = shuffledFacies(unit1);
            shuffledFacies(unit1) = shuffledFacies(unit2);
            shuffledFacies(unit2) = temp;

            %Swap the thicknesses
            temp = shuffledThick(unit1);
            shuffledThick(unit1) = shuffledThick(unit2);
            shuffledThick(unit2) = temp;

            j = j + 1;
        end
        
        infiniteStopCount = infiniteStopCount + 1;
    end
end
