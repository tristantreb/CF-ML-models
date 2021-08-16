function [amIntrDatacube] = amEMMCCreateIntrDatacubeRecovery(amDatacube, amInterventions, measures, align_wind, offset, ninterventions, nmeasures, curveaveragingmethod, datasmoothmethod)

% createIntrDatacube - creates the data cube for offset + alignement window
% by intervention (for each measure)

amIntrDatacube = NaN(ninterventions, align_wind + offset.span-1 + abs(offset.down), nmeasures);
% note: offset.down is present to pour additional data to the right when
% curve is shifted to the left by abs(offset.down) slots

midx = measures.Index(ismember(measures.DisplayName, 'LungFunction'));

for i = 1:ninterventions
    scid   = amInterventions.SmartCareID(i);
    start = amInterventions.IVScaledDateNum(i);
    
    icperiodend = align_wind + offset.span-1 + abs(offset.down); % offset of  means curve cannot be shifted
    %dcperiodend = start-1 + align_wind + offset.up + abs(offset.down);
    dcperiodend = start-1 + align_wind + offset.span-1 + abs(offset.down);
    
    if curveaveragingmethod == 1
        fprintf('*** WARNING *** curveaveragingmethod 1 not implemented');
    else
        icperiodstart = 1; % put 1 because it is used as an index below (should be 0)
        %dcperiodstart = start + offset.down; 
        dcperiodstart = start; 
    end
    
    % ajust to min offset to first datapoint for this intervention
    if dcperiodstart <= 0
        icperiodstart = icperiodstart - dcperiodstart + 1;
        dcperiodstart = 1;
    end

% TODO % why not using the movmax filter for all spirometry measurements?
    for m = 1:nmeasures
        % for datasmoothmethod 2, smooth FEV1 measures with a 3 day max
        % window, else just use raw data
        if (datasmoothmethod == 2 && m == midx)
            amIntrDatacube(i, (icperiodstart:icperiodend), m) = movmax(amDatacube(scid, dcperiodstart:dcperiodend, m), 2, 'omitnan');
        elseif (datasmoothmethod == 3 && m == midx)
            amIntrDatacube(i, (icperiodstart:icperiodend), m) = movmax(amDatacube(scid, dcperiodstart:dcperiodend, m), 3, 'omitnan');    
        else
            amIntrDatacube(i, (icperiodstart:icperiodend), m) = amDatacube(scid, dcperiodstart:dcperiodend, m);
        end
    end
end
            
end

