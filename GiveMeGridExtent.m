function extent = GiveMeGridExtent(timePoint)

switch timePoint
case 1
    % 'E11pt5',
    extent = [70,75,40];
case 2
    % E13pt5',
    extent = [89,109,69];
case 3
    % 'E15pt5',
    extent = [94,132,65];
case 4
    % 'E18pt5',
    extent = [67,43,40];
case 5
    % P4',
    extent = [77,43,50];
case 6
    % 'P14',
    extent = [68,40,50];
case 7
    % 'P28',
    extent = [73,41,53];
end

end
