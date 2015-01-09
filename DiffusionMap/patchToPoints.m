

function points = patchToPoints(patch)

points = squeeze(reshape(patch, size(patch,1)*size(patch,2), 1, 3));

