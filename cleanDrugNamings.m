function out = cleanDrugNamings(str)
% inputs: array, output: array with cleaned drugtherapy names
    out = replace(str,["symkevi","Symkevi Modulator"],"Symkevi");
    out = replace(out,["Kaftrio","Trikaftor","Triple Therapy",...
        "Triple therapy","kaftrio" ],"Trikafta");
    out = replace(out,"ivacaftor","Ivacaftor");
    fprintf('Cleaned drug therapy namings:\n')
    disp(unique(out))
end
