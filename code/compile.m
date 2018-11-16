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






