function s = fixedlengthint(num,le,align)

% generate string s that contains the integer num
% and always has total length le.
% if align = 'l', will left align num
% if align = 'r', will right align num

s  = num2str(num,'%i');
ls = length(s);
if ls > le
    fprintf('warning: number longer than le. Printing full number.')
    return;
end

if ~ischar(align)
    align = 'l';
end

sp = spcs(le-ls);

if strcmp(align,'l');
    s = [s,sp];
else
    s = [sp,s];
end