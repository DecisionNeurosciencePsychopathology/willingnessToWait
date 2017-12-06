function [rects] = setRects(origin)
% sets stimulus positions for qTask_v3

% first is the position of the reward indicator (the "light")
rects.light.rect = [origin(1)-50, origin(2)-300, origin(1)+50, origin(2)-200];

% identify the center of the light for placing text
rects.light.center = [mean(rects.light.rect([1,3])), mean(rects.light.rect([2,4]))];

% next the rect for the left-hand button ("smaller-sooner")
rects.optSmall.rect = [origin(1)-203, origin(2)-100, origin(1)-3, origin(2)+100];

% rect for the right-hand button ("larger-later")
rects.optLarge.rect = [origin(1)+3, origin(2)-100, origin(1)+203, origin(2)+100];

% identify the center of each button (for placing text)
rects.optSmall.center = [mean(rects.optSmall.rect([1,3])), mean(rects.optSmall.rect([2,4]))];
rects.optLarge.center = [mean(rects.optLarge.rect([1,3])), mean(rects.optLarge.rect([2,4]))];

% thick borders for when each is chosen
rects.optSmall.border = rects.optSmall.rect + 6*[-1 -1 1 1];
rects.optLarge.border = rects.optLarge.rect + 6*[-1 -1 1 1];

% text at the bottom of the screen
rects.time = [origin(1)-100, origin(2)+150];
rects.earnings = [origin(1)-100, origin(2)+190];
