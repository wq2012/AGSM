if is_octave()
    mex force_field/GVF.cpp -o force_field/GVF;
    mex force_field/BoundMirrorEnsure.cpp -o force_field/BoundMirrorEnsure;
    mex force_field/BoundMirrorExpand.cpp -o force_field/BoundMirrorExpand;
    mex force_field/BoundMirrorShrink.cpp -o force_field/BoundMirrorShrink;
    mex force_field/ContourDirectionForce.cpp -o force_field/ContourDirectionForce;
    mex force_field/ContourForceArray.cpp -o force_field/ContourForceArray;
    mex force_field/InitialCircle.cpp -o force_field/InitialCircle;

    mex line_fitting/LineInImage.cpp -o line_fitting/LineInImage;
    mex line_fitting/LineNormalForce.cpp -o line_fitting/LineNormalForce;
    mex line_fitting/LineTorque.cpp -o line_fitting/LineTorque;

    mex circle_fitting/CircleNormalForce.cpp -o circle_fitting/CircleNormalForce;
    mex circle_fitting/CircleHoughTransform.cpp -o circle_fitting/CircleHoughTransform;

    mex arbitrary_ellipse_fitting/EllipseTorque.cpp -o arbitrary_ellipse_fitting/EllipseTorque;

    mex spline_contour_fitting/SplineContourForce.cpp -o spline_contour_fitting/SplineContourForce;
else
    mex force_field/GVF.cpp -outdir force_field;
    mex force_field/BoundMirrorEnsure.cpp -outdir force_field;
    mex force_field/BoundMirrorExpand.cpp -outdir force_field;
    mex force_field/BoundMirrorShrink.cpp -outdir force_field;
    mex force_field/ContourDirectionForce.cpp -outdir force_field;
    mex force_field/ContourForceArray.cpp -outdir force_field;
    mex force_field/InitialCircle.cpp -outdir force_field;

    mex line_fitting/LineInImage.cpp -outdir line_fitting;
    mex line_fitting/LineNormalForce.cpp -outdir line_fitting;
    mex line_fitting/LineTorque.cpp -outdir line_fitting;

    mex circle_fitting/CircleNormalForce.cpp -outdir circle_fitting;
    mex circle_fitting/CircleHoughTransform.cpp -outdir circle_fitting;

    mex arbitrary_ellipse_fitting/EllipseTorque.cpp -outdir arbitrary_ellipse_fitting;

    mex spline_contour_fitting/SplineContourForce.cpp -outdir spline_contour_fitting;
end
