function defaultarg(varargin)

[hasval, flagind] = ismember('-value', varargin);

if hasval
    default_val = evalin('caller', varargin{flagind+1});
    varargin = varargin(setdiff(1:end, flagind+[0,1]));
else
    default_val = [];
end

for i=1:numel(varargin)
    if ~evalin('caller', ['exist(''' varargin{i} ''', ''var'')'])
        assignin('caller', varargin{i}, default_val);
    end
end