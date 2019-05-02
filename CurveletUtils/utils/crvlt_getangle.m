function ang = crvlt_getangle(idx, totnr)

% crvlt_getangle(idx, totnr)
%
% This function gives the angles for oscillating direction of curvelets 
% with direction indices idx
% when the total number of directions on this level is totnr

idxinq = mod(idx-1,totnr/4);     % index in quadrant
quadrant = floor((idx-1) / (totnr/4));

slope = ((totnr/8) - idxinq - 0.5)/(totnr/8);

ang = atan(slope) + pi/2 - pi/2 * quadrant;

