function h=plott3(varargin) 
% A 'plot3' wrapper for 3xN array
%Y.M.

if nargin>=3 && isnumeric(varargin{2}) && isnumeric(varargin{3}) && ...
        all(size(varargin{1})==size(varargin{2})) && all(size(varargin{1})==size(varargin{3}))
    h = feval(@plot3, varargin{:});
else    
    xyz = squeeze(varargin{1});
    if isempty(xyz), return; end
    if size(xyz, 2)==3 && size(xyz,1)~=3
        disp('Warning: xyz in plott3 should be a 3*N matrix...');
        xyz=xyz';
    end    
    h = feval(@plot3, xyz(1,:),xyz(2,:),xyz(3,:), varargin{2:end});
end

if nargout==0
    clear h;
end

