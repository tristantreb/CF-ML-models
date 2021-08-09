function [amInterventions] = RamEMMCCalcConfidenceBounds(overall_pdoffset, amInterventions, offset_span, ninterventions, confidencethreshold, confidencemode)

% amEMMCCalcConfidenceBounds - calculates the lower and upper confidence bounds
% around the predicted offset

% note - maximum probability mode only works for prob distributions with a
% max of two peaks.

intrcount = array2table((1:ninterventions)');
intrcount.Properties.VariableNames{1} = 'IntrNbr';
amInterventions = [intrcount, amInterventions];
amInterventions.LowerBound1(:) = amInterventions.Offset;
amInterventions.UpperBound1(:) = amInterventions.Offset;
%amInterventions.LowerBound2(:) = -1;
%amInterventions.UpperBound2(:) = -1;
amInterventions.ConfidenceProb(:) = 0;

for i = 1:ninterventions
    
    adjconfthreshold = confidencethreshold * sum(overall_pdoffset(amInterventions.LatentCurve(i), i, 1:offset_span));
    
    % Contiguous mode: populate the confidence bounds until they contain 
    % at least the probability defined by the confidence threshold. Start 
    % with the most likely day and increment by one day at a time (either 
    % moving left or right to the next most likely
    if confidencemode == 1
        fprintf('Not Implemented');
    
    elseif confidencemode == 2
        % Maximum mode: select days in order of most likely days
        [sortedprob, sortedidx] = sort(overall_pdoffset(amInterventions.LatentCurve(i), i, 1: offset_span),'descend');
        cumprob = 0;
        n = 0;
        while cumprob < adjconfthreshold && n <= size(sortedidx,3)
            n = n + 1;
            cumprob = cumprob + sortedprob(n);
        end
        confidenceidx = sort(sortedidx(1:n),'ascend');
        confidencedays = offsetval(confidenceidx);
        amInterventions.LowerBound1(i) = confidencedays(1);
        amInterventions.UpperBound1(i) = confidencedays(end);
        
%         if amInterventions.LowerBound2(i) == -1
%             amInterventions.UpperBound1(i) = confidencedays(a) - 1;
%         else
%             amInterventions.UpperBound2(i) = confidencedays(a) - 1;
%         end
        amInterventions.ConfidenceProb(i) = cumprob;
        fprintf('Intervention %2d, Latent Curve %1d, Offset = %2d, Lower1 = %2d, Upper1 = %2d, cumprob %.4f\n', ...
            i, amInterventions.LatentCurve(i), amInterventions.Offset(i), amInterventions.LowerBound1(i), amInterventions.UpperBound1(i), ...
            cumprob);%amInterventions.LowerBound2(i), amInterventions.UpperBound2(i), cumprob);
    else
        fprintf('Invalid Confidence mode\n');
    end

end

end
