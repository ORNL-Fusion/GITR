if perDiffusionCoefficientInterpolator_number ==0
    perDiffusionCoefficientInterpolator = @gitrInterpScalar0D;
else if perDiffusionCoefficientInterpolator_number == 1
        perDiffusionCoefficientInterpolator = @gitrInterpScalar1D;
    else if perDiffusionCoefficientInterpolator_number == 2
            perDiffusionCoefficientInterpolator = @gitrPerDiffusionAnalytic;
        else if perDiffusionCoefficientInterpolator_number == 3
                perDiffusionCoefficientInterpolator = @gitrInterpScalar3D;
            end
        end
    end
end


if perDiffusionCoefficientInterpolator_number ==0
    perDiffusionCoeff(:) = perDiffusionCoeff_in;
else if perDiffusionCoefficientInterpolator_number ==1
                perDiffusionCoeff = reshape(dlmread('perDiffCoeff_out.txt','\t'),nXv, nYv, nZv);
    else if perDiffusionCoefficientInterpolator_number ==2
            error('no analytical model for perpendicular diffusion implemented')
        else if perDiffusionCoefficientInterpolator_number ==3
                perDiffusionCoeff = reshape(dlmread('perDiffCoeff_out.txt','\t'),nXv, nYv, nZv);
            end
        end
    end
end


interpolators = {EfieldInterpolator; BfieldInterpolator;...
    FlowVelocityInterpolator; temperatureInterpolator; ...
    densityInterpolator; perDiffusionCoefficientInterpolator};