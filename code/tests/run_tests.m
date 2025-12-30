function run_tests()
    fprintf('Starting AGSM Unit Test Suite...\n');
    
    try
        test_math();
        test_force_field();
        test_line();
        test_circle();
        test_ellipse();
        test_spline();
        
        fprintf('\n============================\n');
        fprintf('ALL TESTS PASSED SUCCESSFULLY\n');
        fprintf('============================\n');
        exit(0);
    catch err
        fprintf('\n============================\n');
        fprintf('TEST SUITE FAILED\n');
        fprintf('Error: %s\n', err.message);
        fprintf('============================\n');
        exit(1);
    end
end
