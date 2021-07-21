function z = div_gfq(q,varargin)
global div_result 

z =varargin{1};

for i=2:length(varargin)
    if varargin{i}>=q
        ME = MException('GfValErr:illegalval', ...
            'illegal value for gf(q)',varargin{i});
        throw(ME)
    end
    z =div_result(z+1,varargin{i}+1);
end

