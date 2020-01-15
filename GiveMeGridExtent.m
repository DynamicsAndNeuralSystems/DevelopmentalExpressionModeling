function extent = GiveMeGridExtent(timePoint)
% Can get it from the main repository as:
coOrds = getCoOrds('wholeBrain','P28');
max(coOrds)-min(coOrds);

switch timePoint
case 1
    % 'E11pt5' (resolution 0.08)
    extent = [44,59,15]*0.08;
case 2
    % E13pt5',
    extent = [57,33,12]*0.1;
case 3
    % 'E15pt5',
    extent = [55,31,14]*0.12;
case 4
    % 'E18pt5',
    extent = [56,29,15]*0.14;
case 5
    % P4',
    extent = [66,34,19]*0.16;
case 6
    % 'P14',
    extent = [64,33,21]*0.2;
case 7
    % 'P28',
    extent = [68,37,21]*0.2;
end

end
