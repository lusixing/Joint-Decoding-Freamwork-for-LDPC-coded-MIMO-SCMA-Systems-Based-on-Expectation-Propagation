function z = add_gfq(q,varargin)
global add_result
z =0;

for i=1:length(varargin)
    
    if length(z)>1 && length(varargin{i})>1
        for j=1:length(varargin{i})
            if varargin{i}(j)>=q
                ME = MException('GfValErr:illegalval', ...
                    'illegal value for gf(q)',varargin{i}(j));
                throw(ME)
            end
            z(j)=add_result(z(j)+1,varargin{i}(j)+1);
        end
    else
        if varargin{i}>=q
            ME = MException('GfValErr:illegalval', ...
                'illegal value for gf(q)',varargin{i});
            throw(ME)
        end
        z =add_result(z+1,varargin{i}+1);
    end
end

