#!/bin/bash -l

#distributed computing toolbox does not seem to be in 2015b
#module load matlab/r2015a    THIS FAILS WITH LICENSE ISSUES!
module load matlab/r2014b

## !! execute on /cfs/klemming/..  !!

nohup matlab -nosplash -nodesktop << EOF > submit_log.txt &
configCluster

% Configure job: email, memory, partition, Project name...
ClusterInfo.setEmailAddress('hufnagel@kth.se')
% request fat node
ClusterInfo.setMemUsage('2000000')
% set allocation
ClusterInfo.setProjectName('2015-16-46')
ClusterInfo.setRequireExclusiveNode(true)
ClusterInfo.setProcsPerNode(48)
ClusterInfo.setWallTime('03:00:00')
ClusterInfo.setUserDefinedOptions(' --output=logfile.%J.out --error=errfile.%J.err -J matlab_pod')

c=parcluster;

j = c.batch('pod');

% Wait for the job to finish before fetching results
j.wait 


%c.getDebugLog(j.Tasks(1)) is sometimes helpful helpful..

% Now that the job has completed, fetch the results
if strcmp(j.State,'finished')
  diary(j) 
  j.fetchOutputs{:}
end

% delete the job
j.delete

disp('end');
EOF
