% Divide an image, get rid of INF

function img = imdiv(im1)

warning off;

img = 1./im1;

img(isinf(img))=0;
img(isnan(img))=0;

warning on;