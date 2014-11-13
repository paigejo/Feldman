function [memStats] = linMem()
if ~isunix
    memStats.error = true;
    return
end

[sts,msg] = unix('free -m | grep Mem:');

if sts %error
    memStats.error = true;
    return
else %no error
    mems = cell2mat(textscan(msg,'%*s %u %u %u %*u %*u %*u','delimiter',' ','collectoutput',true,'multipleDelimsAsOne',true));
    memStats.freeMB = mems(3);
    memStats.usedMB = mems(2);
    memStats.totalMB = mems(1);
    memStats.error = false;
end

end