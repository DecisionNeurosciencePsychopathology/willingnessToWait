function [rew] = RewFunctionWTW()
%Calculate the expected reward for a policy of quitting at time t.
%From McGuire 2012 Rt = 0.15(pt) + 0.01(1-pt)
%Where:
%   pt = proportion of rewards delivered earlier than t

