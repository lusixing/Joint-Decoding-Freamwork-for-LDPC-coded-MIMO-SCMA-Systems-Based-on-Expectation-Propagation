function z = mul_gfq(q,varargin)

global mul_result
z =1;

for i=1:length(varargin)
    if varargin{i}>=q
        ME = MException('GfValErr:illegalval', ...
            'illegal value for gf(q)',varargin{i});
        throw(ME)
    end
    z =mul_result(z+1,varargin{i}+1);
end





