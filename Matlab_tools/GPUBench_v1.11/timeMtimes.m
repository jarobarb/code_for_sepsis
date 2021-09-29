function timeMtimes()
highOverheadWarningId = 'parallel:gpu:gputimeit:HighOverhead';
cleanup = onCleanup( @() warning( 'on', highOverheadWarningId ) );
warning( 'off', highOverheadWarningId );

N = 32;
data = rand(N,N,'gpuArray');
gputimeit( @() data*data )

for ii=1:10
    iTimeit( 'GPU', @() data*data )
end

%-------------------------------------------------------------------------%
function t = iTimeit( device, f )
% Depending by what device is selected, time f with gputimeit(gpu) or
% timeit(cpu)
if strcmp( device, 'GPU' )
    t = gputimeit( f );
else
    t = timeit( f );
end
