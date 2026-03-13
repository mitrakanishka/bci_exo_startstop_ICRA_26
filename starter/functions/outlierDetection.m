function logicalIdx = outlierDetection(data, p)

madfactor = -1 /(sqrt(2)*erfcinv(3/2));  %~1.4826
center = median(data,1,'omitnan');
amad = madfactor*median(abs(data - center), 1, 'omitnan');

lowerbound = center - p*amad;
upperbound = center + p*amad;

logicalIdx = (data > upperbound | data < lowerbound);
