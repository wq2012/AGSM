function assert_equal(actual, expected, tol, name)
    if nargin < 3 || isempty(tol)
        tol = 1e-6;
    end
    if nargin < 4
        name = 'unnamed test';
    end
    
    if all(size(actual) == size(expected)) && all(abs(actual(:) - expected(:)) < tol)
        fprintf('PASS: %s\n', name);
    else
        fprintf('FAIL: %s\n', name);
        fprintf('  Actual: %s\n', mat2str(actual));
        fprintf('  Expected: %s\n', mat2str(expected));
        error('Assertion failed for %s', name);
    end
end
