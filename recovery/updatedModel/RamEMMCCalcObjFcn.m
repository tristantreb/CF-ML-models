function [dist, count, hstg, vshift, isOutlier] = RamEMMCCalcObjFcn(meancurvemean, meancurvestd, amIntrCube, ...
    amHeldBackcube, vshift, isOutlier, outprior, measuresmask, measuresrange, normstd, hstg, currinter, ...
    curroffsetidx, offset, align_wind, nmeasures, update_histogram, sigmamethod, smoothingmethod, allowvshift, vshiftmax)

% amEMMCCalcObjFcn - calculates residual sum of squares distance for points in
% curve vs a given meancurve incorporating offset

dist = 0;
count = 0;
tempmean = zeros(align_wind + offset.span-1, nmeasures);
tempstd  = zeros(align_wind + offset.span-1, nmeasures);
offsetval = offset.down:offset.up;
curroffsetval = offsetval(curroffsetidx); % convert offsetidx to offsetval

if (update_histogram == 1)
    for m = 1:nmeasures
        hstg(1, m, currinter, curroffsetidx) = 0;
    end
end

% check right size
for m = 1:nmeasures
    if smoothingmethod == 2
        tempmean(:,m) = smoothdata(meancurvemean(1, :, m),5);
        tempstd(:,m) = smoothdata(meancurvestd(1, :, m),5);
    else
        tempmean(:,m) = meancurvemean(1, :, m);
        tempstd(:,m) = meancurvestd(1, :, m);
    end
end

% partially vectorised implementation for improved performance
for m = 1:nmeasures
    outdist = -log(outprior) + log(measuresrange(m));
    actidx = ~isnan(amIntrCube(currinter, 1 : align_wind, m)) & amHeldBackcube(currinter, 1 : align_wind, m) == 0;
    offsettempmean = tempmean(abs(offset.down) + 1 + curroffsetval : abs(offset.down) + align_wind + curroffsetval, m);
    intrdata       = amIntrCube(currinter, 1 : align_wind, m);
    
    if allowvshift && sum(actidx) ~= 0
        vshift(1, currinter, m, curroffsetidx) = (sum(offsettempmean(actidx)) - sum(intrdata(actidx))) / sum(actidx);
        if vshift(1, currinter, m, curroffsetidx) > vshiftmax
            vshift(1, currinter, m, curroffsetidx) = vshiftmax;
        elseif vshift(1, currinter, m, curroffsetidx) < -vshiftmax
            vshift(1, currinter, m, curroffsetidx) = -vshiftmax;
        end
    end
    
    if sigmamethod == 4
        regdist = ((offsettempmean - intrdata' - vshift(1, currinter, m, curroffsetidx)) .^ 2) ...
                    ./ (2 * (tempstd(abs(offset.down) + 1 + curroffsetval : abs(offset.down) + align_wind + curroffsetval, m) .^ 2)) ...
                    + log(tempstd(abs(offset.down) + 1 + curroffsetval : abs(offset.down) + align_wind + curroffsetval, m)) ...
                    + log((2 * pi) ^ 0.5) ...
                    - log(1 - outprior);
    else
        regdist = ( (offsettempmean - intrdata' - vshift(1, currinter, m, curroffsetidx)) .^ 2) ...
                    / ( 2 * normstd(currinter, m) ) ...
                    + log(normstd(currinter, m)) ...
                    + log((2 * pi) ^ 0.5) ...
                    - log(1 - outprior);
    end
    
%    % partially vectorised implementation for improved performance
% for m = 1:nmeasures
%     outdist = -log(outprior) + log(measuresrange(m));
%     actidx = ~isnan(amIntrCube(currinter, abs(offset.down) + 1 : abs(offset.down) + align_wind, m)) & amHeldBackcube(currinter, abs(offset.down) + 1 : abs(offset.down) + align_wind, m) == 0;
%     offsettempmean = tempmean(abs(offset.down) + 1 + curroffsetval : abs(offset.down) + align_wind + curroffsetval, m);
%     intrdata       = amIntrCube(currinter, abs(offset.down) + 1 :  abs(offset.down) + align_wind, m);
%     
%     if allowvshift && sum(actidx) ~= 0
%         vshift(1, currinter, m, curroffsetidx) = (sum(offsettempmean(actidx)) - sum(intrdata(actidx))) / sum(actidx);
%         if vshift(1, currinter, m, curroffsetidx) > vshiftmax
%             vshift(1, currinter, m, curroffsetidx) = vshiftmax;
%         elseif vshift(1, currinter, m, curroffsetidx) < -vshiftmax
%             vshift(1, currinter, m, curroffsetidx) = -vshiftmax;
%         end
%     end
%     
%     if sigmamethod == 4
%         regdist = ((offsettempmean - intrdata' - vshift(1, currinter, m, curroffsetidx)) .^ 2) ...
%                     ./ (2 * (tempstd(abs(offset.down) + 1 + curroffsetval : abs(offset.down) + align_wind + curroffsetval, m) .^ 2)) ...
%                     + log(tempstd(abs(offset.down) + 1 + curroffsetval : abs(offset.down) + align_wind + curroffsetval, m)) ...
%                     + log((2 * pi) ^ 0.5) ...
%                     - log(1 - outprior);
%     else
%         regdist = ( (offsettempmean - intrdata' - vshift(1, currinter, m, curroffsetidx)) .^ 2) ...
%                     / ( 2 * normstd(currinter, m) ) ...
%                     + log(normstd(currinter, m)) ...
%                     + log((2 * pi) ^ 0.5) ...
%                     - log(1 - outprior);
%     end

    regdist = regdist(actidx);
    outdist = outdist * ones(sum(actidx), 1);
    ptidx = regdist < outdist;
    totdist = sum(regdist(ptidx)) + sum(outdist(~ptidx));

    if measuresmask(m) == 1
        dist = dist + totdist;
        count = count + sum(actidx);
    end

    if (update_histogram == 1)
        hstg(1, m, currinter, curroffsetidx) = totdist;
    end

end

%fprintf('Intr %2d, Offset %2d: Dist: Old %8.4f, New %8.4f   Count: Old %3d New %3d\n', currinter, curroffset, dist, dist2, count, count2);

end
