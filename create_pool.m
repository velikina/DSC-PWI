function pool=create_pool(nw, isIncr)

%  pool=create_pool(nw, isIncr)
%
% Closes existing pool if its number of workers not equal nw,
% then creates pool with specified number of workers nw. If the number of
% specified workes exceed the maximum number, then creates the pool with
% the maximum possible number of workers. If isIncr > 0, changes the number
% of workers only if it is higher than in the existing pool
% 
% 
% 

if ~check_var('isIncr')
    isIncr=0;
end

pool=gcp('nocreate');

defaultProfile = parallel.defaultClusterProfile();
myCluster = parcluster(defaultProfile);

if nw>myCluster.NumWorkers
    nw=myCluster.NumWorkers;
    if nw>20
        nw=20;
    end
    disp(['Changing the supplied number of workers to max available (' num2str(nw)  ')']);
end

if isempty(pool)
    pool=parpool(nw);
elseif pool.NumWorkers~=nw && ~(isIncr && pool.NumWorkers>nw)        
    delete(gcp);
    pool=parpool(nw);
else
    display('The pool satisfying the inputs is already active');
end

end

