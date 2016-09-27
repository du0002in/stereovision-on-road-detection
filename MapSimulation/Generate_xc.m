load 'map_cordinate4_dotted';
load state.mat % variable S
% rand('state',S);
rng(S);

nSegmentsLD=length(map_x)/300 %300 data points corresponding to about 15meters
segmentLengthLD=round(400+rand(1,nSegmentsLD)*(500-300));
csH = cumsum(segmentLengthLD);
nSegmentsLD = length(find(cumsum(segmentLengthLD)<length(map_x))) ;

% extend length of last segment to complete the road
segmentLengthLD(nSegmentsLD) = length(map_x) - sum(segmentLengthLD(1:nSegmentsLD-1)); 
segmentLengthLD = segmentLengthLD(1:nSegmentsLD);
segmentLD=zeros(1,nSegmentsLD,);
segmentLD(1) = (simulated_lane_width/2) ; % begin at the center of the left lane
segmentLD(2:end)=simulated_lane_width/2+maximumLateralDisplacement*2*(rand(1,nSegmentsLD-1)-0.5);



%-------------------------------------
% Computation of the camera trajectory
%-------------------------------------
% maximum possible number of segments of constant lateral displacement
nSegmentsLD = roadLength/minSegmentLengthLD ; 
segmentLengthLD = round(minSegmentLengthLD + ...
                  rand(1,nSegmentsLD)*(maxSegmentLengthLD - minSegmentLengthLD));
csH = cumsum(segmentLengthLD);
nSegmentsLD = length(find(cumsum(segmentLengthLD)<roadLength)) ;
% extend length of last segment to complete the road
segmentLengthLD(nSegmentsLD) = roadLength - sum(segmentLengthLD(1:nSegmentsLD-1)); 
segmentLengthLD = segmentLengthLD(1:nSegmentsLD);
segmentLD(1) = (laneWidth/2) ; % begin at the center of the right lane
for i=2:nSegmentsLD
    % keep always in the lane to the right of the road center (two lanes=>right lane)
    segmentLD(i) = (laneWidth/2) + maximumLateralDisplacement*2*(rand(1)-0.5) ; 
end
ld = [] ;
for i=1:nSegmentsLD
    ld  = [ld ; segmentLD(i)*ones(segmentLengthLD(i),1)];
end
% heavy smoothing
t = double(-windowSizeLD/2:windowSizeLD/2);
filter = (1/(sqrt(2*pi)*sigmaLD))*exp(-(t.*t)/(2*sigmaLD*sigmaLD));
ld = smoothVector(ld, filter);