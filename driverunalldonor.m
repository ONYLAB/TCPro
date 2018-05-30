for traj = 72:250
    
    [~,~,~,~,Inco(:,traj),ELIS(:,traj),COMBO(:,traj),corrs(:,traj)] = runalldonor(traj);
    
%     load([num2str(traj) '_matlabRcutoff100Prol0597.mat'],'Inco')
%     res(:,traj) = Inco;
end