function enhanced = globalStretching(A, gamma, Tol, alpha)
% Global Stretching
chkCast = class(A);

% Gamma correction
A = A .^ gamma;

% Normalization to contrast stretch to [0,1]
A = mat2gray(A);

% Find limits to stretch the image
clipLimit = stretchlim(A, Tol);

% Adjust the cliplimits
clipLimit = clipLimit + alpha * (max(clipLimit, mean(clipLimit, 2)) - clipLimit);

% Adjust the image intensity values to new cliplimits
enhanced = imadjust(A, clipLimit);
enhanced = cast(enhanced, chkCast);

end