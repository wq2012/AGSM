% Copyright (C) 2012 Quan Wang <wangq10@rpi.edu>
% Standard library check for Octave environment.

function r = is_octave()
% IS_OCTAVE Check if the current environment is GNU Octave.
%    r = is_octave() returns true if the current environment is Octave, 
%    false otherwise.

    persistent x;
    if isempty(x)
        x = (exist('OCTAVE_VERSION', 'builtin') > 0);
    end
    r = x;
end
