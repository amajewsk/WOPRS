function hhmmss_particle_new=fix_SPEC_times(infile,ncfile,diodesize,processflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function fixes SPEC times in clock cycles before converting clock
% cycles to HHMMSS data that can be compared to buffer times. Diode size in
% mm (0.01 for 2DS). Process flag is either debug to view statistics of the
% timing differences and comparison plots, but not write out the corrected 
% information *or* 'write' to run in process mode with no plots and
% overwrite the 'Time_in_seconds' variable with the SPEC particle time and
% update the variable attributes to reflect.
%
% Ex. Usage:    s_compare = fix_SPEC_times('~/Documents/WOPRS_test/SNOWIE/PROC.170119b.2DS.cdf','~/Documents/WOPRS_test/SNOWIE/20170119b.c1.nc',.01,'debug')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Time and buffer data read
    %Load PROC data
    f = netcdf.open(infile,'nowrite');
    buffer_hhmmss = netcdf.getVar(f,netcdf.inqVarID(f,'Time'));
    buffer_msec = netcdf.getVar(f,netcdf.inqVarID(f,'msec'));
    buffer_hhmmss_fine = buffer_hhmmss+buffer_msec/1d3;
    s_particle = netcdf.getVar(f,netcdf.inqVarID(f,'Time_in_seconds'));
    parent_rec_num = netcdf.getVar(f,netcdf.inqVarID(f,'parent_rec_num'));
    clocks = netcdf.getVar(f,netcdf.inqVarID(f,'particle_microsec'));
    slices = netcdf.getVar(f,netcdf.inqVarID(f,'slicecount'));
    netcdf.close(f);
    % %Load HK data
    % parts = strsplit(infile,'/');
    % cdir = strjoin(parts(1:end-1),'/');
    % nameparts = strsplit(parts{end},'.');
    % hkfile = [cdir '/DIMG.' nameparts{2} '.2DS.HK.cdf'];
    % f = netcdf.open(hkfile);
    % hk_hh = netcdf.getVar(f,netcdf.inqVarID(f,'hour','double'));
    % hk_mm = netcdf.getVar(f,netcdf.inqVarID(f,'minute'));
    % hk_ss = netcdf.getVar(f,netcdf.inqVarID(f,'second'));
    % hk_ms = netcdf.getVar(f,netcdf.inqVarID(f,'millisec'));
    % hk_hhmmss = 1d4*double(hk_hh)+1d2*double(hk_mm)+double(hk_ss)+double(hk_ms)*1d-3;
    %Load AC TAS data
    nc = netcdf.open(ncfile,'nowrite');
    nc_time = netcdf.getVar(f,netcdf.inqVarID(nc,'time'));
    nc_hhmmss = sec2hhmmss(mod(nc_time,86400));
    try
        tas = netcdf.getVar(f,netcdf.inqVarID(nc,'tas'),'double');
    catch 
        tas = netcdf.getVar(f,netcdf.inqVarID(nc,'TAS'),'double');
    end
    netcdf.close(nc);

    % Rolls first
    mask_roll = find(clocks(2:end)-clocks(1:end-1) < -2^31); % (tw2 - img_last(2)) < 0 && (tw1 - img_last(1)) < 0 && abs(tw1 - HK_last(1)) > 2^15 && HK_last(1) < 60000
    for i = 1:numel(mask_roll)
        clocks(mask_roll(i):end) = clocks(mask_roll(i):end) + 2^32 - 1;
    end
    
    %Recursively fix negative delta's
    clocks = negative_times(clocks);

    %Propogate times forward from buffer starts using TAS
    rec_rolls = find(diff(parent_rec_num)<0);
    for k = 1:numel(rec_rolls)
        %fix buffer number rolls
        if k < numel(rec_rolls)
            parent_rec_num(rec_rolls(k)+1:rec_rolls(k+1)) = parent_rec_num(rec_rolls(k))+(parent_rec_num(rec_rolls(k)+1:rec_rolls(k+1))); 
        else
            parent_rec_num(rec_rolls(k)+1:end) = parent_rec_num(rec_rolls(k))+(parent_rec_num(rec_rolls(k)+1:end)); 
        end
    end
    [recs,i_recs,~] = unique(parent_rec_num);
    delta_buffers = zeros(numel(recs),1);
    cum_adjs = delta_buffers;
    hhmmss_particle_new = zeros(size(s_particle));
    ia_clocks = zeros(numel(s_particle),1);
    temp_buffer_time = buffer_hhmmss_fine(end);
    last_j = numel(s_particle);
    for i = numel(i_recs):-1:2% loop over all unique buffer record numbers from last to first
        i_tas = nc_hhmmss == buffer_hhmmss(i_recs(i));% index of matching tas time to current buffer
        tas_temp = max([50.0 tas(i_tas)]);
        cum_adj = 0;
        delta_buffer = hhmmss2sec(temp_buffer_time) - hhmmss2sec(buffer_hhmmss_fine(i_recs(i)-1));
        for j = last_j:-1:i_recs(i)
            %s_off = hhmmss2sec(buffer_hhmmss_fine(j))-s_particle(j);
            if j == last_j % Last (in time) in buffer, first fixed
                delta = max([slices(j) 0]);
            else
                delta = clocks(j)-clocks(j-1)+max([slices(j) 0]);% delta in clock cycles from last TW (no delta for first buffer)
            end
            adj = max([delta/tas_temp*(diodesize/(10^3)) 0]); %convert clocks to seconds (= clock_delta / TAS * mks_diode_resolution)
            cum_adj = cum_adj + adj;
            if j == last_j && delta_buffer > 0
                new_hhmmss = buffer_hhmmss_fine(i_recs(i));% last particle in buffer and buffer delta > 0
            else
                new_hhmmss = sec2hhmmss(hhmmss2sec(hhmmss_particle_new(j+1))-adj);%max([(hhmmss_particle_new(j+1) - adj) buffer_hhmmss_fine(i_recs(i)-1)]);
            end
            %disp(['Delta of ',num2str(adj),', cumulatively ',num2str(cum_adj)]);% 
            hhmmss_particle_new(j) = new_hhmmss;
            ia_clocks(j) = adj;% if saving inter arrival times in seconds;
        end
        % Catch any PbP rollbacks after reconstructed times (current buffer)
        ind_back = i_recs(i) - 1 + find(hhmmss_particle_new(i_recs(i):last_j) < buffer_hhmmss_fine(i_recs(i)-1));
        if buffer_hhmmss_fine(i_recs(i)-1) <= buffer_hhmmss_fine(i_recs(i)) %if buffer times in order
            hhmmss_particle_new(ind_back) = buffer_hhmmss_fine(i_recs(i)-1);
        else %fix rollbacks on last reconstructed particle time in order
            hhmmss_particle_new(ind_back) = hhmmss_particle_new(ind_back(end)+1);
        end
        disp(['Buffer ' num2str(parent_rec_num(last_j-1)) ' end ',num2str(buffer_hhmmss_fine(i_recs(i))),', start ',num2str(buffer_hhmmss_fine(i_recs(i)-1)),', total delta of ',num2str(delta_buffer)]);
        disp(['Accumulated particle deltas ',num2str(100.0*cum_adj/delta_buffer),'% of ',num2str(last_j-i_recs(i)+1),' buffer delta']);
        temp_buffer_time = buffer_hhmmss_fine(i_recs(i)-1);
        last_j = i_recs(i)-1;
        delta_buffers(i) = delta_buffer;
        cum_adjs(i) = cum_adj;
    end
    error = delta_buffers - cum_adjs;
    mse = immse(delta_buffers,cum_adjs);

    % Time Conversions
    s_buffer_fine=hhmmss2sec(buffer_hhmmss_fine);
    s_particle_new=hhmmss2sec(hhmmss_particle_new);
    s_offset = s_buffer_fine - s_particle(i_recs(2)); %offset between particle seconds (s from start of flight)
    s_particle_day = s_particle + s_offset; %...converted to HHMMSS (s from start of day)
    
    if strcmp(processflag, 'debug') 
        %% Diagnostic Plots
        disp(['The MSE of the buffer timing words appears to be ',num2str(mse),' s']);
        [N_b,bin_b] = histcounts(delta_buffers,linspace(0,1.,1001));
        [N_p,bin_p] = histcounts(cum_adjs,linspace(0,1.,1001));
        tiledlayout(2,1)
        % Time delta (buffer/PbP) vs. Log10 count freq
        ax1 = nexttile;
        plot(bin_p(2:end),log10(N_p),'g',bin_b(2:end),log10(N_b),'r')
        xlabel(ax1,'Elapsed Time (s, upper bin edge)')
        ylabel(ax1,'Log_{10}(Counts)')
        % Record Number vs. Time delta (buffer/PbP)
        ax2 = nexttile;
        plot(parent_rec_num,s_particle_new,'-+g',parent_rec_num,s_buffer_fine,'r') %parent_rec_num,s_particle_day,'-+b',
        xlabel(ax2,'Buffer Number')
        ylabel(ax2,'Seconds from Start of Day')
        ytickformat(ax2,'%6.2f')
        
        grid on
        grid minor
        hold off
    end

    if strcmpi(processflag, 'WRITE')
        f = netcdf.open(infile,upper(processflag));
        v_ptime = netcdf.inqVarID(f,'Time_in_seconds');
        v_ia = netcdf.inqVarID(f,'inter_arrival');
        netcdf.putAtt(f,v_ptime,'long_name','SPEC reconstructed particle time')
        netcdf.putAtt(f,v_ptime,'units','HHMMSS.SSS (UTC)')
        netcdf.putAtt(f,v_ia,'units','seconds')
        netcdf.putVar(f,v_ptime,hhmmss_particle_new);
        netcdf.putVar(f,v_ia,ia_clocks);
        v_global = netcdf.getConstant('NC_GLOBAL');
        hist = netcdf.getAtt(f,v_global,'history');
        netcdf.reDef(f);
        netcdf.putAtt(f,v_global,'history',[hist,' ',char(datetime('now')),': particle times reconstructed']);
        netcdf.close(f);
    end
end

function corrected = negative_times(clocks)
    % Fix TW jumps (singular)
    ind_neg = find(diff(clocks) < 0); % indices of bad times (not first, not last)
    diff_about = clocks(ind_neg+1)-clocks(ind_neg-1);
    single_up = find((5d3 > diff_about) & (diff_about > 0));
    clocks(ind_neg(single_up)) = .5*(clocks(ind_neg(single_up)-1)+clocks(ind_neg(single_up)+1));
    ind_neg = find(diff(clocks) < 0);
    %plot(ind_neg,clocks(ind_neg),'r+')

    %tic
    % Fix longer tw jumps
    for i = 1:numel(ind_neg)
        i_n = ind_neg(i)+1;% index of the n-th (last bad point)
        i_0 = find(clocks < clocks(i_n));% index of the last good point
        i_0 = i_0(end);

        % Fix any TW jumps (after unrolling, TW decreases discretely but not for only 1 mistimed cluster)
        % Check to see if jump is once a
        if abs(.25*pi - atan((clocks(i_0+1)-clocks(i_0))/(clocks(i_n-1)-clocks(i_n)))) > pi/8
            %offset all prior clocks to align
            clocks(i_n:end) = clocks(i_n:end)+clocks(i_n-1)-clocks(i_n);% offset remaining clocks forward to next good clock
        else
            n_bad = i_n-i_0-1;% num elevated timing words
            clock_diff = mean([-clocks(i_n)+clocks(i_n-1) clocks(i_0+1)-clocks(i_0)]);
            clocks(i_0+1:i_n-1)  = clocks(i_0+1:i_n-1) - clock_diff;
            %b = clocks(i_0) + ((clocks(i_n)-clocks(i_0))/(n_bad))*linspace(1,n_bad,n_bad);%*linspace(i_0+1,i_n-1,n_bad);
        end
    end
    %toc

    ind_neg = find(diff(clocks) < 0);
    if numel(ind_neg) > 0
        clocks = negative_times(clocks);
        disp('recursed');
    end
    corrected = clocks;
end

