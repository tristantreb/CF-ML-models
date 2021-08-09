function [amIntrDatacube] = amEMMCCreateIntrDatacubeRecovery_withoutoffset(amDatacube, amInterventions, measures, align_wind, max_offset, ninterventions, nmeasures, curveaveragingmethod, datasmoothmethod)

% createIntrDatacube - creates the data cube for offset + alignement window
% by intervention (for each measure)

amIntrDatacube = NaN(ninterventions, max_offset + align_wind - 1, nmeasures);
midx = measures.Index(ismember(measures.DisplayName, 'LungFunction'));

for i = 1:ninterventions
    scid   = amInterventions.SmartCareID(i);
    start = amInterventions.IVScaledDateNum(i);
    
    icperiodend = align_wind + (max_offset - 1); % max_offset of 1 means curve cannot be shifted
    dcperiodend = start + align_wind + (max_offset - 1) - 1; % -1 to have the same amount of days
    
    if curveaveragingmethod == 1
        fprintf('*** WARNING *** curveaveragingmethod 1 not implemented');
    else
        icperiodstart = 1; % put 1 because it is used as an index below (should be 0)
        dcperiodstart = start;
    end
    
%     if dcperiodstart <= 0 % TODO % what happens if it is greater than the
%     last measurement recorded
%         icperiodstart = icperiodstart - dcperiodstart + 1;
%         dcperiodstart = 1;
%     end
    
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

