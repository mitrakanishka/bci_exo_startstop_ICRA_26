function [data_filt, zf_out] = func_pseudo_online_filter(data, param, zf_in, Type)
% non-causal filter for filtering a data frame (e.g., 32 data points)
% inputs:
%  - param: matlab struct containing parameters for the desired filtering type
%     - Type == HPF: [param.hpa, param.hpb] = butter(...,'high')
%     - Type == LPF: [param.lpa, param.lpb] = butter(...,'low')
%     - Type == AVG: param.avg is the window size (in seconds) for moving average smoothing
%     - Type == EMA: param.ema is the alpha parameter of EMA: mu{t} = alpha*x{t} + (1-alpha)*mu{t-1}. alpha in (0,1)
%     - Type == CAR: set param to []
%  - data: data to be filtered. If multichannel, organize channels as columns (i.e., size = Nt * Nch)
%  - zf_in: final filter conditions from previous frame. If first frame, provide zf_in = [] !! 
%  - Type:
%     - HPF: high-pass filtering
%     - LPF: low-pass filtering
%     - AVG: moving average filtering (low-pass)
%     - EMA: exponential moving average filtering (high-pass effect, removes the EMA from the signal!!)
%     - CAR: common average reference spatial filter
% 
% outputs:
%  - data_filt: filtered data
%  - zf_out: final filter conditions
% 
% Last modified: 12/07/2022 by F. Samuel Racz {fsr324@austin.utexas.edu}

if strcmp(Type, 'HPF')
    if isempty(zf_in)
        [data_filt, zf_out] = filter(param.hpa, param.hpb, data);
    else
        [data_filt, zf_out] = filter(param.hpa, param.hpb, data, zf_in);
    end
elseif strcmp(Type, 'LPF')
    if isempty(zf_in)
        [data_filt, zf_out] = filter(param.lpa, param.lpb, data);
    else
        [data_filt, zf_out] = filter(param.lpa, param.lpb, data, zf_in);
    end
elseif strcmp(Type, 'BPF')
    if isempty(zf_in)
        [data_filt, zf_out] = filter(param.bpa, param.bpb, data);
    else
        [data_filt, zf_out] = filter(param.bpa, param.bpb, data, zf_in);
    end
elseif strcmp(Type, 'AVG')
    if isempty(zf_in)
        [data_filt, zf_out] = filter(ones(1, ceil(param.avg*param.fs))./param.avg/param.fs, 1, data);
    else
        [data_filt, zf_out] = filter(ones(1, ceil(param.avg*param.fs))./param.avg/param.fs, 1, data, zf_in);
    end
elseif strcmp(Type, 'EMA')
    data_filt = zeros(size(data));
    for t = 1:size(data,1)
        if isempty(zf_in)
            zf_out = data(1,:);
        else
            zf_out = param.ema*data(t,:) + (1-param.ema)*zf_in;
        end
        data_filt(t,:) = data(t,:) - zf_out;
        zf_in = zf_out;
    end
elseif strcmp(Type, 'CAR')
    data_filt = data - mean(data,2);
    zf_out = [];
end


end