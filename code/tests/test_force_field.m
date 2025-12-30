% Copyright (C) 2012 Quan Wang <wangq10@rpi.edu>,
% Signal Analysis and Machine Perception Laboratory,
% Department of Electrical, Computer, and Systems Engineering,
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA
%
% You are free to use this software for academic purposes if you cite our paper:
% Quan Wang, Kim L. Boyer,
% The active geometric shape model: A new robust deformable shape model and its applications,
% Computer Vision and Image Understanding, Volume 116, Issue 12, December 2012, Pages 1178-1194,
% ISSN 1077-3142, 10.1016/j.cviu.2012.08.004.
%
% For commercial use, please contact the authors.

function test_force_field()
    fprintf('\nTesting Force Field Functions...\n');
    addpath('..');
    addpath('../force_field');
    
    % create small test image
    I = zeros(20, 20);
    I(10, 10) = 100;
    I = gaussianBlur(I, 2);
    
    % test GVF
    [u, v] = GVF(I, 1, 0.1, 10);
    assert_equal(size(u), [20, 20], [], 'GVF size u');
    assert_equal(size(v), [20, 20], [], 'GVF size v');
    
    % test InitialCircle
    [xc, yc, r] = InitialCircle(I);
    assert_equal(isscalar(xc) && isscalar(yc) && isscalar(r), true, [], 'InitialCircle output types');
    
    % test ContourForceArray
    x = [10 11 12];
    y = [10 10 10];
    f = ContourForceArray(u, v, x, y);
    assert_equal(size(f), [3, 2], [], 'ContourForceArray size');

    % test BoundMirrorEnsure
    I_padded = BoundMirrorEnsure(I);
    assert_equal(size(I_padded), size(I), [], 'BoundMirrorEnsure size');

    % test BoundMirrorExpand
    I_expanded = BoundMirrorExpand(I);
    assert_equal(size(I_expanded), size(I) + 2, [], 'BoundMirrorExpand size');

    % test BoundMirrorShrink
    I_shrunk = BoundMirrorShrink(I_expanded);
    assert_equal(size(I_shrunk), size(I), [], 'BoundMirrorShrink size');

    % test ContourDirectionForce
    dx = ones(20, 20);
    dy = zeros(20, 20);
    p = [1 0];
    dir_f = ContourDirectionForce(dx, dy, x, y, p);
    assert_equal(size(dir_f), [1, 1], [], 'ContourDirectionForce size');
    
    fprintf('Force field tests completed successfully.\n');
end
