% MULTI_WAITBAR Multimodal waitbar
%h = multi_waitbar(X, 'message') behaves the same way as the regular waitbar
% function is the global variable g_waitbar_type is empty or set to
% 'graphical'. Otherwise, it will silently do nothing.

function h = multi_waitbar(varargin)
    global g_waitbar_type;

    if isempty(g_waitbar_type) || strcmp(g_waitbar_type, 'graphical')
        if isinf(varargin{1})
            close(varargin{2});
        else
            h = waitbar(varargin{:});
        end
    else
        h = [];
    end
end

