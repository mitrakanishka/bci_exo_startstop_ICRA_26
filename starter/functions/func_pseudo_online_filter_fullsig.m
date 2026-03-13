function [data_filt, zf_out] = func_pseudo_online_filter_fullsig(param,data,zf_in,Type)

data_filt = zeros(size(data));

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