function measures = getMeasuresMask(measuresmask, measures)
 % reset measures mask
 measures.Mask(:) = 0;
 % set the measures mask depending on option chosen
 if measuresmask == 0
     % all
     measures.Mask(:) = 1;
 elseif measuresmask(:) == 1
     % all except activity
     idx = ~ismember(measures.DisplayName, {'Activity'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 2
     % cough, lung function and wellness
     idx = ismember(measures.DisplayName, {'Cough', 'LungFunction', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 3
     % all except activity and lung function
     idx = ~ismember(measures.DisplayName, {'Activity', 'LungFunction'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 4
     % all except activity and weight
     idx = ~ismember(measures.DisplayName, {'Activity', 'Weight'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 5
     % all except weight
     idx = ~ismember(measures.DisplayName, {'Weight'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 6
     % project breathe
     idx = ismember(measures.DisplayName, {'Calorie', 'Cough', 'FEV1', 'MinsAsleep', 'MinsAwake', ...
         'O2Saturation', 'PulseRate', 'RestingHR', 'Temperature', 'Weight', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 7
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'FEV1', 'MinsAwake', 'O2Saturation', ...
         'RestingHR', 'Temperature', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 8
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'FEF2575', 'MinsAsleep', 'O2Saturation', ...
         'RestingHR', 'Temperature', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 9
     % project climb
     idx = ismember(measures.DisplayName, {'Appetite', 'Cough', 'SleepActivity', ...
         'SputumVolume', 'Tiredness', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 10
     % project climb
     idx = ~ismember(measures.DisplayName, {'Activity', 'Temperature'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 11
     % project climb
     idx = ismember(measures.DisplayName, {'Appetite', 'Cough', 'O2Saturation', 'PulseRate', ...
         'SleepActivity', 'SputumVolume', 'Tiredness', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 12
     % project climb
     idx = ismember(measures.DisplayName, {'Appetite', 'Cough', 'FEV1', 'O2Saturation', 'PulseRate', ...
         'SleepActivity', 'SputumVolume', 'Tiredness', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 13
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'FEV1', 'MinsAsleep', ...
         'RestingHR', 'Temperature', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 14
     % project climb
     idx = ismember(measures.DisplayName, {'Appetite', 'Cough', 'FEV1', 'O2Saturation', 'PulseRate', ...
         'SleepActivity', 'SputumVolume', 'Temperature', 'Tiredness', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 15
     % project climb
     idx = ismember(measures.DisplayName, {'Appetite', 'Cough', 'FEV1', 'O2Saturation', 'PulseRate', ...
         'SleepActivity', 'SputumVolume', 'Temperature', 'Tiredness', 'Weight', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 16
     % project climb
     idx = ismember(measures.DisplayName, {'SputumVolume'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 17
     % project climb
     idx = ismember(measures.DisplayName, {'Appetite', 'Cough', 'O2Saturation', 'PulseRate', ...
         'SleepActivity', 'SputumVolume', 'Temperature', 'Tiredness', 'Weight', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 18
     % project climb
     idx = ismember(measures.DisplayName, {'Appetite', 'Cough', 'InterpFEV1', 'O2Saturation', 'PulseRate', ...
         'SleepActivity', 'SputumVolume', 'Temperature', 'Tiredness', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 19
     % project climb
     idx = ismember(measures.DisplayName, {'Appetite', 'Cough', 'InterpFEV1', 'O2Saturation', 'PulseRate', ...
         'SleepActivity', 'SputumVolume', 'Temperature', 'Tiredness', 'InterpWeight', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 20
     % project climb
     idx = ismember(measures.DisplayName, {'Appetite', 'Cough', 'O2Saturation', 'PulseRate', ...
         'SleepActivity', 'SputumVolume', 'Temperature', 'Tiredness', 'InterpWeight', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 21
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'FEV1', 'MinsAsleep', 'O2Saturation', ...
         'RestingHR', 'Temperature', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 22
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'FEV1', 'MinsAsleep', 'O2Saturation', ...
         'PulseRate', 'Temperature', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 23
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'FEV1', 'MinsAsleep', 'O2Saturation', ...
         'RestingHR', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 24
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'FEV1', 'MinsAsleep', 'MinsAwake', 'O2Saturation', ...
         'RestingHR', 'Wellness'});
     measures.Mask(idx) = 1;    
 elseif measuresmask == 25
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'FEV1', 'MinsAsleep', 'MinsAwake', 'O2Saturation', ...
         'RestingHR', 'Temperature', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 26
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'LungFunction', 'MinsAsleep', ...
         'PulseRate', 'Temperature', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 27
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'LungFunction', 'MinsAsleep', 'O2Saturation', ...
         'PulseRate', 'Temperature', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 28
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'LungFunction', 'MinsAsleep', ...
         'RestingHR', 'Temperature', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 29
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'LungFunction', 'MinsAsleep', 'O2Saturation', ...
         'RestingHR', 'Temperature', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 30
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'FEV1', 'FEF2575', 'O2Saturation', ...
         'RestingHR', 'PulseRate', 'Weight', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 31
     % project breathe
     idx = ismember(measures.DisplayName, {'FEV1'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 32
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'FEV1', 'FEF2575', 'O2Saturation', ...
         'PulseRate', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 33
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough', 'MinsAsleep', 'FEV1', 'FEF2575', 'O2Saturation', ...
         'PulseRate', 'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 34
     % project breathe - reference for recovery
     idx = ismember(measures.DisplayName, {'Cough', 'FEV1', 'FEF2575', 'MinsAsleep', 'O2Saturation', ...
         'PulseRate', 'Wellness', 'Temperature'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 35
     % project breathe
     idx = ismember(measures.DisplayName, {'Cough'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 36
     % project breathe
     idx = ismember(measures.DisplayName, {'Wellness'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 37
     % project breathe
     idx = ismember(measures.DisplayName, {'MinsAsleep'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 38
     % project breathe
     idx = ismember(measures.DisplayName, {'O2Saturation'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 39
     % project breathe
     idx = ismember(measures.DisplayName, {'PulseRate'});
     measures.Mask(idx) = 1;
 elseif measuresmask == 40
     % project breathe
     idx = ismember(measures.DisplayName, {'Temperature'});
     measures.Mask(idx) = 1;
 else
     % shouldn't ever get here - but default to just cough if it ever
     % happens
     idx = ismember(measures.DisplayName, {'Cough'});
     measures.Mask(idx) = 1;
end
end
